/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2026-08-31 $

// Massimo Petracca - ASDEA Software, Italy
//
// J2 plasticity with Armstrong-Frederick / Chaboche kinematic hardening and
// IMPL-EX. See ASDSteel3DMaterial.h for the model, the proof that its uniaxial
// restriction is exactly ASDSteel1DMaterial, and the Voigt convention.

#include "ASDSteel3DMaterial.h"
#include "PlaneStrainMaterial.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <classTags.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// ==================================================================== //
//  1  Voigt utilities                                                  //
// ==================================================================== //

namespace {

	// EVERY 6-VECTOR IN THIS FILE CARRIES THE TENSOR COMPONENTS OF THE SHEAR.
	// The engineering convention lives in setTrialStrain, getStrain, the two
	// tangents and the plastic-strain response, and nowhere else. The one
	// consequence to keep in mind is that a full contraction doubles the last
	// three terms, which is what contract() is for.
	//
	// These five are the same utilities ASDPlasticDamageConcrete3DMaterial keeps
	// in its own anonymous namespace. They are duplicated rather than shared: a
	// header with them in it would be worth having, but factoring them out of a
	// file that works today is a change to that file, and it is not this one.

	const double SQ23 = std::sqrt(2.0 / 3.0);
	const double SQ32 = std::sqrt(1.5);

	enum ErrorCodes {
		EC_Generic = -1,
		EC_IMPLEX_Error_Control = -10
	};

	inline double traceOf(const Vector& a) { return a(0) + a(1) + a(2); }

	// a : b, the full contraction of two symmetric tensors
	inline double contract(const Vector& a, const Vector& b) {
		return a(0) * b(0) + a(1) * b(1) + a(2) * b(2)
			+ 2.0 * (a(3) * b(3) + a(4) * b(4) + a(5) * b(5));
	}

	inline void deviatorOf(const Vector& a, Vector& out) {
		double p = traceOf(a) / 3.0;
		out(0) = a(0) - p;
		out(1) = a(1) - p;
		out(2) = a(2) - p;
		out(3) = a(3);
		out(4) = a(4);
		out(5) = a(5);
	}

	// sigma = lam*tr(e)*I + 2*mu*e, in closed form, on TENSOR strain. mu2 = 2*mu.
	// This is the dividend of the tensor-component convention: the elastic law
	// needs no 6x6 matrix and no engineering factor in either direction.
	inline void elasticStress(double lam, double mu2, const Vector& e, Vector& out) {
		double lam2 = lam + mu2;
		out(0) = lam2 * e(0) + lam * e(1) + lam * e(2);
		out(1) = lam * e(0) + lam2 * e(1) + lam * e(2);
		out(2) = lam * e(0) + lam * e(1) + lam2 * e(2);
		out(3) = mu2 * e(3);
		out(4) = mu2 * e(4);
		out(5) = mu2 * e(5);
	}

	// The 6x6 elastic operator, TENSOR in / TENSOR out: its shear columns carry
	// the factor two of the contraction, which is what makes a plain matrix
	// product compose correctly. getInitialTangent()/getTangent() are the two
	// places it is brought back to the engineering-input convention.
	inline void elasticMatrix(double lam, double mu2, Matrix& out) {
		out.Zero();
		double lam2 = lam + mu2;
		out(0, 0) = out(1, 1) = out(2, 2) = lam2;
		out(0, 1) = out(1, 0) = out(0, 2) = out(2, 0) = out(1, 2) = out(2, 1) = lam;
		out(3, 3) = out(4, 4) = out(5, 5) = mu2;
	}

	// TO THE ENGINEERING-INPUT CONVENTION: halve the shear COLUMNS. A 6x6 built
	// here acts on tensor components and OpenSees hands it engineering ones. Call
	// it ONCE, at the very end - it is not the same thing as the factor two a
	// contraction carries, and applying both to the same object cancels them.
	inline void toEngineeringInput(Matrix& D) {
		for (int i = 0; i < 6; ++i)
			for (int j = 3; j < 6; ++j)
				D(i, j) *= 0.5;
	}

	// (1 - exp(-x))/x, continued to 1 at x = 0.
	//
	// THE POINT OF WRITING IT THIS WAY. The exact integral of the AF law over a
	// step is a_i = A_i*n + (a_i_commit - A_i*n)*exp(-gamma_i*Dp) with
	// A_i = sqrt(2/3)*H_i/gamma_i, which divides by gamma_i and therefore has to
	// branch when a backstress does not saturate. Written through this function
	// the same expression is Phi_i = sqrt(2/3)*H_i*Dp*g(gamma_i*Dp), which at
	// gamma_i = 0 gives sqrt(2/3)*H_i*Dp - linear Prager hardening, the correct
	// limit, reached continuously and with no special case anywhere.
	inline double gfun(double x) {
		if (std::abs(x) < 1.0e-4)
			return 1.0 - x * (0.5 - x * (1.0 / 6.0 - x / 24.0));
		return -std::expm1(-x) / x;
	}

	bool string_to_double(const std::string& text, double& num) {
		num = 0.0;
		try {
			num = std::stod(text);
			return true;
		}
		catch (...) {
			return false;
		}
	}

} // namespace

// ==================================================================== //
//  2  the parser                                                       //
// ==================================================================== //

void* OPS_ASDSteel3DMaterial(void)
{
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDSteel3D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	static const char* msg =
		"nDMaterial ASDSteel3D $tag $E $nu $sy $su $eu "
		"<-rho $rho> "
		"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> "
		"<-implexAbort> <-implexAlpha $alpha> "
		"<-tangent> <-fullTangent> <-elasticTangent>";

	if (OPS_GetNumRemainingInputArgs() < 6) {
		opserr << "nDMaterial ASDSteel3D Error: few arguments (< 6).\n" << msg << "\n";
		return nullptr;
	}

	int numData = 1;
	int tag = 0;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "nDMaterial ASDSteel3D Error: invalid 'tag'.\n" << msg << "\n";
		return nullptr;
	}

	double E = 0.0, nu = 0.0, sy = 0.0, su = 0.0, eu = 0.0, rho = 0.0;
	auto lam_get = [&numData](double* val, const char* name) -> bool {
		if (OPS_GetDouble(&numData, val) != 0) {
			opserr << "nDMaterial ASDSteel3D Error: invalid '" << name << "'.\n" << msg << "\n";
			return false;
		}
		return true;
		};
	if (!lam_get(&E, "E")) return nullptr;
	if (!lam_get(&nu, "nu")) return nullptr;
	if (!lam_get(&sy, "sy")) return nullptr;
	if (!lam_get(&su, "su")) return nullptr;
	if (!lam_get(&eu, "eu")) return nullptr;

	auto lam_optional_double = [&numData](const char* name, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() < 1) {
			opserr << "nDMaterial ASDSteel3D Error: '" << name << "' requested but not provided.\n" << msg << "\n";
			return false;
		}
		if (OPS_GetDouble(&numData, &value) < 0) {
			opserr << "nDMaterial ASDSteel3D Error: failed to get '" << name << "'.\n" << msg << "\n";
			return false;
		}
		return true;
		};

	ASDSteel3DMaterial::InputParameters params;

	// THE PARSING LOOP HAS A CATCH-ALL, and that is the one difference from the
	// 1D's loop worth pointing at: that one is a chain of independent ifs with no
	// final else, so a misspelled option is swallowed in silence and the user
	// gets the default behaviour under the impression they asked for something
	// else.
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-rho") == 0) {
			if (!lam_optional_double("rho", rho)) return nullptr;
		}
		else if (strcmp(value, "-implex") == 0) {
			params.implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			params.implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "nDMaterial ASDSteel3D Error: '-implexControl' given without the next 2 "
					"arguments $implexErrorTolerance $implexTimeReductionLimit.\n" << msg << "\n";
				return nullptr;
			}
			if (!lam_optional_double("implexErrorTolerance", params.implex_error_tolerance)) return nullptr;
			if (!lam_optional_double("implexTimeReductionLimit", params.implex_time_redution_limit)) return nullptr;
		}
		else if (strcmp(value, "-implexAbort") == 0) {
			params.implex_abort_on_error = true;
		}
		else if (strcmp(value, "-implexAlpha") == 0) {
			if (!lam_optional_double("implexAlpha", params.implex_alpha)) return nullptr;
		}
		else if (strcmp(value, "-tangent") == 0) {
			params.tangent_type = ASDSteel3DMaterial::Tangent_Numerical;
		}
		else if (strcmp(value, "-fullTangent") == 0) {
			params.tangent_type = ASDSteel3DMaterial::Tangent_Full;
		}
		else if (strcmp(value, "-elasticTangent") == 0) {
			params.tangent_type = ASDSteel3DMaterial::Tangent_Elastic;
		}
		else {
			// a number here means a positional argument too many, which is a
			// different mistake from a misspelled keyword and deserves to be
			// named as such
			double dummy;
			if (string_to_double(value, dummy))
				opserr << "nDMaterial ASDSteel3D Error: unexpected value '" << value
				<< "'. This material takes 5 values after the tag ($E $nu $sy $su $eu).\n" << msg << "\n";
			else
				opserr << "nDMaterial ASDSteel3D Error: unknown option '" << value << "'.\n" << msg << "\n";
			return nullptr;
		}
	}

	// checks. Two of these the 1D does NOT do and they are not oversights to
	// inherit: it validates E only against zero, and it divides by eu with no
	// guard at all, so eu = 0 silently produces an infinite H1. ASDConcrete3D
	// for its part never validates nu
	if (E <= 0.0) {
		opserr << "nDMaterial ASDSteel3D Error: invalid 'E' (" << E << "). It should be strictly positive.\n";
		return nullptr;
	}
	if (nu <= -1.0 || nu >= 0.5) {
		opserr << "nDMaterial ASDSteel3D Error: invalid 'nu' (" << nu << "). It should be in (-1, 0.5); "
			"at 0.5 the bulk modulus is infinite.\n";
		return nullptr;
	}
	if (sy <= 0.0) {
		opserr << "nDMaterial ASDSteel3D Error: invalid 'sy' (" << sy << "). It should be strictly positive.\n";
		return nullptr;
	}
	if (su <= sy) {
		opserr << "nDMaterial ASDSteel3D Error: invalid 'su' (" << su << "). It should be larger than 'sy' ("
			<< sy << ").\n";
		return nullptr;
	}
	if (eu <= sy / E) {
		opserr << "nDMaterial ASDSteel3D Error: invalid 'eu' (" << eu << "). It should be larger than the yield "
			"strain sy/E (" << sy / E << ").\n";
		return nullptr;
	}
	if (rho < 0.0) {
		opserr << "nDMaterial ASDSteel3D Error: invalid 'rho' (" << rho << "). It should be non-negative.\n";
		return nullptr;
	}
	if (params.implex_error_tolerance <= 0.0 || params.implex_time_redution_limit <= 0.0) {
		opserr << "nDMaterial ASDSteel3D Error: the '-implexControl' values should be strictly positive.\n";
		return nullptr;
	}

	// CHABOCHE CALIBRATION FROM (E, sy, su, eu). Copied from
	// ASDSteel1DMaterial.cpp:2178-2200 and NOT re-derived: the proof in the
	// header shows H_i and gamma_i map one-to-one between the 1D and the 3D with
	// no rescaling, so any drift here would break the equivalence that is this
	// material's main verification.
	//
	// Worth stating a fact the 1D leaves implicit: H1/gamma1 = 0.9*(su - sy) and
	// H2/gamma2 = 0.1*(su - sy) IDENTICALLY, independent of n and of eu. So
	// sy + H1/g1 + H2/g2 == su exactly by construction, which is why
	// stressReference() can recover su from the parameters instead of storing
	// it; eu and n set only the saturation RATE, and alpha and m the split
	// between the two time scales.
	{
		const double n = 400.0;
		const double m = 50.0;
		const double alpha = 0.9;
		double dy_norm = (su - sy) / E;
		double H1_norm = (1.0 / n) / eu;
		double gamma1 = H1_norm / dy_norm;
		double H1 = H1_norm * E;
		double H2 = H1 * m;
		double gamma2 = gamma1 * m;
		params.H1 = H1 * alpha;
		params.gamma1 = gamma1;
		params.H2 = H2 * (1.0 - alpha);
		params.gamma2 = gamma2;
	}
	params.E = E;
	params.nu = nu;
	params.rho = rho;
	params.sy = sy;

	NDMaterial* instance = new ASDSteel3DMaterial(tag, params);
	if (instance == nullptr) {
		opserr << "nDMaterial ASDSteel3D Error: failed to allocate a new material.\n";
		return nullptr;
	}
	return instance;
}

// ==================================================================== //
//  3  life-cycle                                                       //
// ==================================================================== //

ASDSteel3DMaterial::ASDSteel3DMaterial(int _tag, const ASDSteel3DMaterial::InputParameters& _params)
	: NDMaterial(_tag, ND_TAG_ASDSteel3DMaterial)
	, params(_params)
{
	C = getInitialTangent();
}

ASDSteel3DMaterial::ASDSteel3DMaterial()
	: NDMaterial(0, ND_TAG_ASDSteel3DMaterial)
{
}

ASDSteel3DMaterial::~ASDSteel3DMaterial()
{
}

double ASDSteel3DMaterial::getRho(void)
{
	return params.rho;
}

// ==================================================================== //
//  4  the implicit return mapping                                      //
// ==================================================================== //

int ASDSteel3DMaterial::integrate(void)
{
	const double mu2 = params.E / (1.0 + params.nu);
	const double lam = params.nu * mu2 / (1.0 - 2.0 * params.nu);
	const double G2 = mu2;                 // 2*G
	const double sy23 = SQ23 * params.sy;

	static Vector de(6);
	static Vector s_tr(6);
	static Vector xi_tr(6);
	static Vector n(6);
	static Vector w(6);
	static Vector wp(6);
	static Vector xi(6);
	static Vector aux(6);

	// EVERY PASS RESTARTS FROM THE COMMITTED STATE. The 1D needs an exception to
	// this on the IMPL-EX correction pass, because its series/RVE Newton has to
	// start from the delivered slip strain; there is no inner Newton here, every
	// pass is a pure function of (committed state, strain), and the rule is the
	// simpler one. It is also what makes the material safe inside
	// PlaneStressMaterial, whose static condensation calls setTrialStrain many
	// times per outer call.
	ep = ep_commit;
	a1 = a1_commit;
	a2 = a2_commit;
	p = p_commit;
	plastic = false;
	failed = false;
	dgamma = 0.0;
	n_flow.Zero();

	// elastic predictor
	for (int i = 0; i < 6; ++i)
		de(i) = strain(i) - ep_commit(i);
	elasticStress(lam, mu2, de, stress);
	double pm = traceOf(stress) / 3.0;
	deviatorOf(stress, s_tr);
	for (int i = 0; i < 6; ++i)
		xi_tr(i) = s_tr(i) - a1_commit(i) - a2_commit(i);

	double r = std::sqrt(std::max(0.0, contract(xi_tr, xi_tr)));
	if (r - sy23 <= 0.0) {
		// elastic. stress already holds the answer
		tg_r = r;
		return 0;
	}

	for (int i = 0; i < 6; ++i)
		n(i) = xi_tr(i) / r;

	// the scalar Newton, on Dp - the accumulated equivalent plastic strain
	// increment - rather than on Dgamma, so the two tolerances are literally the
	// 1D's and the extrapolated quantity is the same one the 1D extrapolates
	constexpr int MAX_ITER = 1000;
	constexpr double F_REL_TOL = 1.0e-6;   // on the residual, relative to sy
	constexpr double L_ABS_TOL = 1.0e-8;   // absolute, on Dp

	double Dp = 0.0;
	double E1 = 1.0, E2 = 1.0, Phi1 = 0.0, Phi2 = 0.0, Psi1 = 0.0, Psi2 = 0.0;
	double beta = r, D = r, Rp = -1.0;
	bool converged = false;

	for (int iter = 0; iter < MAX_ITER; ++iter) {
		double x1 = params.gamma1 * Dp;
		double x2 = params.gamma2 * Dp;
		double g1 = gfun(x1);
		double g2 = gfun(x2);
		Psi1 = x1 * g1;  E1 = 1.0 - Psi1;
		Psi2 = x2 * g2;  E2 = 1.0 - Psi2;
		Phi1 = SQ23 * params.H1 * Dp * g1;
		Phi2 = SQ23 * params.H2 * Dp * g2;

		beta = r - G2 * SQ32 * Dp - (Phi1 + Phi2);
		for (int i = 0; i < 6; ++i)
			w(i) = Psi1 * a1_commit(i) + Psi2 * a2_commit(i);

		// ||xi||, with xi = beta*n + w and n a unit tensor
		double c = contract(n, w);
		D = std::sqrt(std::max(0.0, beta * beta + 2.0 * beta * c + contract(w, w)));
		if (D <= 1.0e-14 * sy23) {
			// the relative stress collapsed: the flow direction is undefined and
			// there is nothing admissible to return
			failed = true;
			break;
		}
		double R = D - sy23;

		for (int i = 0; i < 6; ++i)
			xi(i) = beta * n(i) + w(i);

		// dR/dDp. d(beta)/dDp = -2G*sqrt(3/2) - sqrt(2/3)*sum(H_i*E_i) and
		// d(w)/dDp = sum(gamma_i*E_i*a_i_commit), both exact for the exponential
		// map, so this Newton is the exact one - as the 1D's is, for the same
		// reason
		double bp = -G2 * SQ32 - SQ23 * (params.H1 * E1 + params.H2 * E2);
		for (int i = 0; i < 6; ++i)
			wp(i) = params.gamma1 * E1 * a1_commit(i) + params.gamma2 * E2 * a2_commit(i);
		for (int i = 0; i < 6; ++i)
			aux(i) = (bp * n(i) + wp(i)) / D;
		Rp = contract(xi, aux);

		if (!(Rp < -1.0e-12 * G2)) {
			// no descent: the same guard the 1D writes as 'if (dF == 0) break'
			failed = true;
			break;
		}

		double dDp = -R / Rp;
		Dp += dDp;
		if (Dp < 0.0)
			Dp = 0.0;

		if (std::abs(R) < F_REL_TOL * sy23 && std::abs(dDp) < L_ABS_TOL) {
			converged = true;
			break;
		}
	}

	if (!converged) {
		failed = true;
		return EC_Generic;
	}

	// re-evaluate the frozen-map quantities at the converged Dp and accept
	{
		double x1 = params.gamma1 * Dp;
		double x2 = params.gamma2 * Dp;
		double g1 = gfun(x1);
		double g2 = gfun(x2);
		Psi1 = x1 * g1;
		Psi2 = x2 * g2;
		Phi1 = SQ23 * params.H1 * Dp * g1;
		Phi2 = SQ23 * params.H2 * Dp * g2;
	}
	double Dg = SQ32 * Dp;

	p = p_commit + Dp;
	for (int i = 0; i < 6; ++i) {
		// TENSOR storage: no weight matrix. In engineering storage this same
		// line would need a factor two on the shears
		ep(i) = ep_commit(i) + Dg * n(i);
		a1(i) = a1_commit(i) + Phi1 * n(i) - Psi1 * a1_commit(i);
		a2(i) = a2_commit(i) + Phi2 * n(i) - Psi2 * a2_commit(i);
		// the deviatoric part, assembled so the mean stress is EXACTLY the
		// elastic one: plastic flow here is deviatoric by construction and
		// tr(sigma) must not drift with it
		stress(i) = s_tr(i) - G2 * Dg * n(i);
	}
	stress(0) += pm;
	stress(1) += pm;
	stress(2) += pm;

	dgamma = Dg;
	plastic = true;
	n_flow = n;

	// what computeTangent() needs
	tg_r = r;
	tg_beta = beta;
	tg_Rp = Rp;
	for (int i = 0; i < 6; ++i)
		tg_xih(i) = xi(i) / D;

	return 0;
}

// ==================================================================== //
//  5  the explicit (extrapolated) pass                                 //
// ==================================================================== //

double ASDSteel3DMaterial::timeFactor(void) const
{
	if (dtime_n_commit <= 0.0)
		return 1.0;
	double f = dtime_n / dtime_n_commit * params.implex_alpha;
	if (!(f > 0.0)) {
		// CLAMPED AT ZERO. A negative ratio extrapolates an irreversible process
		// BACKWARDS: the accumulated plastic strain would stop being monotone and
		// the backstresses would run back down their own saturation curve, which
		// is not a state this material can be in. Clamping degrades the step to a
		// purely elastic prediction - wrong by O(dt) like everything else here,
		// but a state that EXISTS
		return 0.0;
	}
	return f;
}

void ASDSteel3DMaterial::extrapolate(void)
{
	const double mu2 = params.E / (1.0 + params.nu);
	const double lam = params.nu * mu2 / (1.0 - 2.0 * params.nu);

	static Vector de(6);

	double Dg = timeFactor() * dgamma_commit;
	if (Dg > 0.0) {
		double Dp = SQ23 * Dg;
		double x1 = params.gamma1 * Dp;
		double x2 = params.gamma2 * Dp;
		double g1 = gfun(x1);
		double g2 = gfun(x2);
		double Psi1 = x1 * g1;
		double Psi2 = x2 * g2;
		double Phi1 = SQ23 * params.H1 * Dp * g1;
		double Phi2 = SQ23 * params.H2 * Dp * g2;

		p = p_commit + Dp;
		for (int i = 0; i < 6; ++i) {
			ep(i) = ep_commit(i) + Dg * n_commit(i);
			a1(i) = a1_commit(i) + Phi1 * n_commit(i) - Psi1 * a1_commit(i);
			a2(i) = a2_commit(i) + Phi2 * n_commit(i) - Psi2 * a2_commit(i);
		}
		n_flow = n_commit;
		plastic = true;
	}
	else {
		// no extrapolation to make: either the committed step was elastic - in
		// which case dgamma_commit is zero AND n_commit is zero, the exact
		// counterpart of the 1D's sg_commit == 0 - or the time factor clamped
		p = p_commit;
		ep = ep_commit;
		a1 = a1_commit;
		a2 = a2_commit;
		n_flow.Zero();
		plastic = false;
	}

	for (int i = 0; i < 6; ++i)
		de(i) = strain(i) - ep(i);
	elasticStress(lam, mu2, de, stress);

	dgamma = Dg;
	failed = false;
}

// ==================================================================== //
//  6  the tangent                                                      //
// ==================================================================== //

void ASDSteel3DMaterial::elasticTangent(Matrix& out) const
{
	double mu2 = params.E / (1.0 + params.nu);
	double lam = params.nu * mu2 / (1.0 - 2.0 * params.nu);
	elasticMatrix(lam, mu2, out);
}

const Matrix& ASDSteel3DMaterial::getInitialTangent(void)
{
	static Matrix D(6, 6);
	elasticTangent(D);
	toEngineeringInput(D);
	return D;
}

const Matrix& ASDSteel3DMaterial::getTangent(void)
{
	return C;
}

int ASDSteel3DMaterial::numericalTangent(void)
{
	// Six forward-perturbed solves plus the unperturbed one, which is also the
	// one that leaves the state behind. The perturbation is a fraction of the
	// YIELD STRAIN: ASDConcrete3D takes its own from the hardening law's strain
	// tolerance, and this material has no hardening law to take it from.
	//
	// This is a verification tool and a way out if the analytical tangent is
	// ever found wrong, not something to run an analysis on: it costs seven
	// return mappings per material point per iteration.
	static Matrix Cnum(6, 6);
	static Vector strain0(6);
	const double PERT = 1.0e-7 * params.sy / params.E;

	strain0 = strain;
	for (int j = 0; j < 6; ++j) {
		strain = strain0;
		strain(j) += PERT;
		int rc = integrate();
		if (rc < 0) { strain = strain0; return rc; }
		for (int i = 0; i < 6; ++i)
			Cnum(i, j) = stress(i);
	}
	strain = strain0;
	int rc = integrate();
	if (rc < 0)
		return rc;
	for (int j = 0; j < 6; ++j)
		for (int i = 0; i < 6; ++i)
			Cnum(i, j) = (Cnum(i, j) - stress(i)) / PERT;

	// built by perturbing the TENSOR strain, so it is tensor-in and needs the
	// same halving every other operator in this file gets
	C = Cnum;
	toEngineeringInput(C);
	return 0;
}

void ASDSteel3DMaterial::computeTangent(void)
{
	static Matrix Ce(6, 6);
	elasticTangent(Ce);

	// UNDER IMPL-EX THIS IS NOT AN APPROXIMATION OF ANYTHING. The extrapolated
	// update writes sigma = C_e:(eps - ep) with ep built from dgamma_commit and
	// n_commit, both frozen data of the step and independent of eps, so sigma is
	// AFFINE in eps and C_e IS the algorithmic tangent, exactly. That is the
	// property the whole scheme exists for: a constant, symmetric,
	// positive-definite stiffness for the global Newton.
	if (params.implex || !plastic || failed || params.tangent_type == Tangent_Elastic) {
		C = Ce;
		toEngineeringInput(C);
		return;
	}
	if (tg_r <= 0.0 || !(tg_Rp < 0.0)) {
		C = Ce;
		toEngineeringInput(C);
		return;
	}

	const double mu2 = params.E / (1.0 + params.nu);
	const double G2 = mu2;
	const double K = params.E / (3.0 * (1.0 - 2.0 * params.nu));

	const Vector& n = n_flow;
	double theta = 1.0 - G2 * dgamma / tg_r;
	double br = tg_beta / tg_r;
	// (xi_hat : n), which is 1 exactly when the AF memory term is collinear with
	// the frozen normal - i.e. on every radial path
	double xn = contract(tg_xih, n);

	// P_s*n and P_s*xi_hat as ROW vectors: the factor two the contraction
	// carries. It is NOT the engineering halving, which is applied once at the
	// end to the whole operator; applying both to the same object cancels them
	static Vector Psn(6);
	static Vector Psx(6);
	for (int i = 0; i < 3; ++i) {
		Psn(i) = n(i);
		Psx(i) = tg_xih(i);
	}
	for (int i = 3; i < 6; ++i) {
		Psn(i) = 2.0 * n(i);
		Psx(i) = 2.0 * tg_xih(i);
	}

	// K*m(x)m + 2G*theta*Idev, tensor in / tensor out
	C.Zero();
	double c_dev = G2 * theta;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			C(i, j) = K - c_dev / 3.0;
	for (int i = 0; i < 3; ++i)
		C(i, i) += c_dev;
	for (int i = 3; i < 6; ++i)
		C(i, i) = c_dev;

	// the plastic corrections, both rank one on n
	double c_nn = G2 * G2 * dgamma / tg_r;                 // 4G^2*Dgamma/r
	double c_pl = G2 * G2 * SQ32 / tg_Rp;                  // 4G^2*sqrt(3/2)/R'
	if (params.tangent_type == Tangent_Full) {
		// the full consistent tangent. NON-SYMMETRIC through the n(x)xi_hat term,
		// and not because the normal is frozen: the Armstrong-Frederick recovery
		// term is non-associative, so every consistent tangent of this model is
		double c1 = c_nn + c_pl * xn * (1.0 - br);
		double c2 = c_pl * br;
		for (int i = 0; i < 6; ++i)
			for (int j = 0; j < 6; ++j)
				C(i, j) += c1 * n(i) * Psn(j) + c2 * n(i) * Psx(j);
	}
	else {
		// Tangent_Collinear, the default: xi_hat := n, so xn = 1 and the two
		// brackets collapse onto one rank-one term. SYMMETRIC once the shear
		// columns are halved, exact on radial paths, O(step^2) otherwise
		double c1 = c_nn + c_pl;
		for (int i = 0; i < 6; ++i)
			for (int j = 0; j < 6; ++j)
				C(i, j) += c1 * n(i) * Psn(j);
	}

	toEngineeringInput(C);
}

// ==================================================================== //
//  7  state handling                                                   //
// ==================================================================== //

int ASDSteel3DMaterial::setTrialStrain(const Vector& v)
{
	int retval = 0;

	// this material point takes part in the current step, so it takes part in
	// the IMPL-EX error aggregate (see IMPLEXManager.h)
	implexTouch();

	// save dT
	if (!dtime_is_user_defined) {
		dtime_n = ops_Dt;
		if (!commit_done) {
			dtime_0 = dtime_n;
			dtime_n_commit = dtime_n;
		}
	}

	// ENGINEERING IN, TENSOR INSIDE: one of the only three places the
	// distinction exists - see the note at the top of this file
	strain(0) = v(0);
	strain(1) = v(1);
	strain(2) = v(2);
	strain(3) = 0.5 * v(3);
	strain(4) = 0.5 * v(4);
	strain(5) = 0.5 * v(5);

	// The numerical tangent builds its own operator as a side effect of the
	// seven solves it runs, so it is the one path that must not go through
	// computeTangent() afterwards
	bool tangent_done = false;
	if (params.implex) {
		// under IMPL-EX the tangent is the elastic operator whatever the user
		// asked for, so -tangent and -fullTangent are simply not reachable here
		extrapolate();
	}
	else if (params.tangent_type == Tangent_Numerical) {
		retval = numericalTangent();
		tangent_done = true;
	}
	else {
		retval = integrate();
	}
	if (retval < 0)
		return retval;

	if (!tangent_done)
		computeTangent();

	// RECORD WHAT THIS STEP DELIVERS, in ONE place covering every branch above.
	// Under IMPL-EX commitState re-solves implicitly and installs that answer
	// over 'stress', so without this the extrapolated stress would be gone by
	// the time any recorder could ask for it
	stress_implex = stress;

	if (params.implex && params.implex_control && retval == 0) {
		// THE LEGACY in-material measurement, read through the NON-DESTRUCTIVE
		// peek and not measured inline, so that the two places that want the
		// error cannot drift. The error control does not need this at all: the
		// convergence test wrapper measures once per step, through
		// computeImplexErrorMetric(), instead of once per iteration
		double err = computeImplexErrorMetric();
		// and only if the user asked for the old behaviour, fail here
		if (params.implex_abort_on_error && !(err <= params.implex_error_tolerance)) {
			if (dtime_n >= params.implex_time_redution_limit * dtime_0)
				retval = EC_IMPLEX_Error_Control;
		}
	}

	return retval;
}

int ASDSteel3DMaterial::setTrialStrain(const Vector& v, const Vector& /*r*/)
{
	return setTrialStrain(v);
}

int ASDSteel3DMaterial::setTrialStrainIncr(const Vector& v)
{
	static Vector aux(6);
	// strain is stored with TENSOR shear, the increment arrives engineering
	aux(0) = strain(0) + v(0);
	aux(1) = strain(1) + v(1);
	aux(2) = strain(2) + v(2);
	aux(3) = 2.0 * strain(3) + v(3);
	aux(4) = 2.0 * strain(4) + v(4);
	aux(5) = 2.0 * strain(5) + v(5);
	return setTrialStrain(aux);
}

int ASDSteel3DMaterial::setTrialStrainIncr(const Vector& v, const Vector& /*r*/)
{
	return setTrialStrainIncr(v);
}

const Vector& ASDSteel3DMaterial::getStrain(void)
{
	static Vector strain_eng(6);
	strain_eng(0) = strain(0);
	strain_eng(1) = strain(1);
	strain_eng(2) = strain(2);
	strain_eng(3) = 2.0 * strain(3);
	strain_eng(4) = 2.0 * strain(4);
	strain_eng(5) = 2.0 * strain(5);
	return strain_eng;
}

const Vector& ASDSteel3DMaterial::getStress(void)
{
	return stress;
}

int ASDSteel3DMaterial::commitState(void)
{
	if (params.implex) {
		// what the extrapolated step delivered is already in stress_implex,
		// recorded by setTrialStrain, and integrate() below does not touch it.
		// THE IMPLICIT SOLUTION IS WHAT GETS COMMITTED: the explicit one carried
		// the step, this one carries the state. So this same call measures the
		// error and does the first half of the commit - there is no second
		// return mapping anywhere
		int retval = integrate();
		if (retval < 0)
			return retval;
		implex_error = implexStressGap(stress_implex, stress);

		// FREEZE what the next explicit pass will use. With the normal frozen at
		// the trial relative stress, the direction the return path TOOK over the
		// step and the flow direction at its end are the same vector, so unlike a
		// curved return path there is nothing to choose here
		if (plastic && dgamma > 0.0) {
			n_commit = n_flow;
			dgamma_commit = dgamma;
		}
		else {
			n_commit.Zero();
			dgamma_commit = 0.0;
		}
	}

	// OUTSIDE the implex block, which is where ASDSteel1D and ASDConcrete3D both
	// put it and where it belongs: 'the dt of the last committed step' is true
	// whatever scheme is running. Only timeFactor() reads it, and only under
	// IMPL-EX, so keeping it inside would be inert for the model - but it is
	// also PUBLISHED, as dTimeCommit of the 'time' response, and on the implicit
	// path it would then read a permanent zero. revertToLastCommit() copies it
	// back into dtime_n too
	dtime_n_commit = dtime_n;

	strain_commit = strain;
	stress_commit = stress;
	ep_commit = ep;
	a1_commit = a1;
	a2_commit = a2;
	p_commit = p;
	commit_done = true;

	return 0;
}

int ASDSteel3DMaterial::revertToLastCommit(void)
{
	strain = strain_commit;
	stress = stress_commit;
	ep = ep_commit;
	a1 = a1_commit;
	a2 = a2_commit;
	p = p_commit;
	dgamma = 0.0;
	n_flow.Zero();
	plastic = false;
	failed = false;
	dtime_n = dtime_n_commit;
	return 0;
}

int ASDSteel3DMaterial::revertToStart(void)
{
	strain.Zero();
	strain_commit.Zero();
	stress.Zero();
	stress_commit.Zero();
	stress_implex.Zero();
	ep.Zero();
	ep_commit.Zero();
	a1.Zero();
	a1_commit.Zero();
	a2.Zero();
	a2_commit.Zero();
	p = 0.0;
	p_commit = 0.0;
	dgamma = 0.0;
	dgamma_commit = 0.0;
	n_flow.Zero();
	n_commit.Zero();
	dtime_n = 0.0;
	dtime_n_commit = 0.0;
	dtime_0 = 0.0;
	dtime_is_user_defined = false;
	commit_done = false;
	implex_error = 0.0;
	plastic = false;
	failed = false;
	tg_r = 0.0;
	tg_beta = 0.0;
	tg_Rp = 0.0;
	tg_xih.Zero();
	C = getInitialTangent();
	return 0;
}

// ==================================================================== //
//  8  copy and info                                                    //
// ==================================================================== //

NDMaterial* ASDSteel3DMaterial::getCopy(void)
{
	return new ASDSteel3DMaterial(*this);
}

NDMaterial* ASDSteel3DMaterial::getCopy(const char* code)
{
	if (strcmp(code, "ThreeDimensional") == 0 || strcmp(code, "3D") == 0)
		return getCopy();
	// PLANE STRAIN IS NOT IN THE BASE CLASS'S LIST, but the wrapper exists and is
	// general over any "ThreeDimensional" material, so there is no reason for
	// this material not to answer for it. Plane stress, plate fibre and beam
	// fibre do come from NDMaterial::getCopy below, for free
	if (strcmp(code, "PlaneStrain") == 0 || strcmp(code, "PlaneStrain2D") == 0) {
		NDMaterial* copy = getCopy();
		NDMaterial* clone = new PlaneStrainMaterial(getTag(), *copy);
		delete copy;
		return clone;
	}
	return NDMaterial::getCopy(code);
}

const char* ASDSteel3DMaterial::getType(void) const
{
	return "ThreeDimensional";
}

int ASDSteel3DMaterial::getOrder(void) const
{
	return 6;
}

void ASDSteel3DMaterial::Print(OPS_Stream& s, int flag)
{
	s << "ASDSteel3D Material, tag: " << this->getTag() << "\n";
	s << "  E = " << params.E << "  nu = " << params.nu << "  rho = " << params.rho << "\n";
	s << "  sy = " << params.sy << "  su (implied) = " << stressReference() << "\n";
	s << "  H1 = " << params.H1 << "  gamma1 = " << params.gamma1
		<< "  H2 = " << params.H2 << "  gamma2 = " << params.gamma2 << "\n";
	s << "  implex = " << (params.implex ? 1 : 0) << "\n";
}

// ==================================================================== //
//  9  serialization                                                    //
// ==================================================================== //

// The fixed part of the double payload. Kept next to the two functions that
// fill and drain it, and CHECKED against the cursor at the end of both: the 1D
// twin carried two counts that had drifted from the field lists (19 against 22,
// then 13 against 16), so the tail of its data was written past the end of an
// undersized Vector, in silence. A hand-maintained count cannot be trusted; a
// hand-maintained count that is asserted can.
//
//   12 six-vectors: strain, strain_commit, stress, stress_commit, stress_implex,
//                   ep, ep_commit, a1, a1_commit, a2, a2_commit, n_commit  = 72
//   p, p_commit, dgamma_commit                                            =  3
//   dtime_n, dtime_n_commit, dtime_0                                      =  3
//   implex_error                                                          =  1
//   E, nu, rho, sy, H1, gamma1, H2, gamma2,
//   implex_error_tolerance, implex_time_redution_limit, implex_alpha      = 11
static const int ASDSteel3D_NDATA_D = 90;

int ASDSteel3DMaterial::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;

	// integers: the tag, the flags, and AS THE LAST ENTRY the size of the double
	// payload, so the receiver never has to recompute an agreement it cannot see
	ID idata(8);
	int ic = 0;
	idata(ic++) = getTag();
	idata(ic++) = static_cast<int>(params.implex);
	idata(ic++) = static_cast<int>(params.implex_control);
	idata(ic++) = static_cast<int>(params.implex_abort_on_error);
	idata(ic++) = params.tangent_type;
	idata(ic++) = static_cast<int>(dtime_is_user_defined);
	idata(ic++) = static_cast<int>(commit_done);
	idata(ic++) = ASDSteel3D_NDATA_D;
	res = theChannel.sendID(getDbTag(), commitTag, idata);
	if (res < 0) {
		opserr << "ASDSteel3DMaterial::sendSelf() - failed to send ID\n";
		return res;
	}

	Vector ddata(ASDSteel3D_NDATA_D);
	int c = 0;
	auto put6 = [&ddata, &c](const Vector& v) { for (int i = 0; i < 6; ++i) ddata(c++) = v(i); };
	put6(strain);
	put6(strain_commit);
	put6(stress);
	put6(stress_commit);
	put6(stress_implex);
	put6(ep);
	put6(ep_commit);
	put6(a1);
	put6(a1_commit);
	put6(a2);
	put6(a2_commit);
	put6(n_commit);
	ddata(c++) = p;
	ddata(c++) = p_commit;
	ddata(c++) = dgamma_commit;
	ddata(c++) = dtime_n;
	ddata(c++) = dtime_n_commit;
	ddata(c++) = dtime_0;
	ddata(c++) = implex_error;
	ddata(c++) = params.E;
	ddata(c++) = params.nu;
	ddata(c++) = params.rho;
	ddata(c++) = params.sy;
	ddata(c++) = params.H1;
	ddata(c++) = params.gamma1;
	ddata(c++) = params.H2;
	ddata(c++) = params.gamma2;
	ddata(c++) = params.implex_error_tolerance;
	ddata(c++) = params.implex_time_redution_limit;
	ddata(c++) = params.implex_alpha;

	if (c != ASDSteel3D_NDATA_D) {
		opserr << "ASDSteel3DMaterial::sendSelf() - the payload size (" << ASDSteel3D_NDATA_D
			<< ") does not match the number of fields written (" << c << "). This is a bug in "
			"this material, not in the input.\n";
		return -1;
	}

	res = theChannel.sendVector(getDbTag(), commitTag, ddata);
	if (res < 0) {
		opserr << "ASDSteel3DMaterial::sendSelf() - failed to send Vector\n";
		return res;
	}
	return 0;
}

int ASDSteel3DMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;

	ID idata(8);
	res = theChannel.recvID(getDbTag(), commitTag, idata);
	if (res < 0) {
		opserr << "ASDSteel3DMaterial::recvSelf() - failed to receive ID\n";
		return res;
	}
	int ic = 0;
	setTag(idata(ic++));
	params.implex = static_cast<bool>(idata(ic++));
	params.implex_control = static_cast<bool>(idata(ic++));
	params.implex_abort_on_error = static_cast<bool>(idata(ic++));
	params.tangent_type = idata(ic++);
	dtime_is_user_defined = static_cast<bool>(idata(ic++));
	commit_done = static_cast<bool>(idata(ic++));
	int nv_dbl = idata(ic++);

	Vector ddata(nv_dbl);
	res = theChannel.recvVector(getDbTag(), commitTag, ddata);
	if (res < 0) {
		opserr << "ASDSteel3DMaterial::recvSelf() - failed to receive Vector\n";
		return res;
	}

	int c = 0;
	auto get6 = [&ddata, &c](Vector& v) { for (int i = 0; i < 6; ++i) v(i) = ddata(c++); };
	get6(strain);
	get6(strain_commit);
	get6(stress);
	get6(stress_commit);
	get6(stress_implex);
	get6(ep);
	get6(ep_commit);
	get6(a1);
	get6(a1_commit);
	get6(a2);
	get6(a2_commit);
	get6(n_commit);
	p = ddata(c++);
	p_commit = ddata(c++);
	dgamma_commit = ddata(c++);
	dtime_n = ddata(c++);
	dtime_n_commit = ddata(c++);
	dtime_0 = ddata(c++);
	implex_error = ddata(c++);
	params.E = ddata(c++);
	params.nu = ddata(c++);
	params.rho = ddata(c++);
	params.sy = ddata(c++);
	params.H1 = ddata(c++);
	params.gamma1 = ddata(c++);
	params.H2 = ddata(c++);
	params.gamma2 = ddata(c++);
	params.implex_error_tolerance = ddata(c++);
	params.implex_time_redution_limit = ddata(c++);
	params.implex_alpha = ddata(c++);

	if (c != nv_dbl) {
		opserr << "ASDSteel3DMaterial::recvSelf() - read " << c << " fields out of a payload of "
			<< nv_dbl << ". The two ends of the channel disagree on the layout.\n";
		return -1;
	}

	// A UNIT VECTOR THAT WENT THROUGH A CHANNEL IS NOT EXACTLY A UNIT VECTOR any
	// more, and this one is multiplied by 2G*Dgamma on every extrapolated step
	double nn = std::sqrt(std::max(0.0, contract(n_commit, n_commit)));
	if (nn > 0.0) {
		for (int i = 0; i < 6; ++i)
			n_commit(i) /= nn;
	}

	// derived, never sent
	dgamma = 0.0;
	n_flow.Zero();
	plastic = false;
	failed = false;
	C = getInitialTangent();

	return 0;
}

// ==================================================================== //
// 10  parameters                                                       //
// ==================================================================== //

int ASDSteel3DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// 1000 - elasticity & mass
	if (strcmp(argv[0], "E") == 0) {
		param.setValue(params.E);
		return param.addObject(1000, this);
	}
	if (strcmp(argv[0], "v") == 0 || strcmp(argv[0], "nu") == 0) {
		param.setValue(params.nu);
		return param.addObject(1001, this);
	}
	if (strcmp(argv[0], "rho") == 0) {
		param.setValue(params.rho);
		return param.addObject(1002, this);
	}
	// 2000 - time. Each of these takes the time step away from ops_Dt, which is
	// how a driver imposes a history whose dt is not the pseudo-time step
	if (strcmp(argv[0], "dTime") == 0) {
		param.setValue(dtime_n);
		return param.addObject(2000, this);
	}
	if (strcmp(argv[0], "dTimeCommit") == 0) {
		param.setValue(dtime_n_commit);
		return param.addObject(2001, this);
	}
	if (strcmp(argv[0], "dTimeInitial") == 0) {
		param.setValue(dtime_0);
		return param.addObject(2002, this);
	}
	return -1;
}

int ASDSteel3DMaterial::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
	case 1000:
		params.E = info.theDouble;
		return 0;
	case 1001:
		params.nu = info.theDouble;
		return 0;
	case 1002:
		params.rho = info.theDouble;
		return 0;
	case 2000:
		dtime_n = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
	case 2001:
		dtime_n_commit = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
	case 2002:
		dtime_0 = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
	default:
		break;
	}
	return -1;
}

// ==================================================================== //
// 11  responses                                                        //
// ==================================================================== //

Response* ASDSteel3DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	auto make_resp = [&output, this](int rid, const Vector& v, const std::vector<std::string>* labels = nullptr) -> MaterialResponse* {
		output.tag("NdMaterialOutput");
		output.attr("matType", getClassType());
		output.attr("matTag", getTag());
		if (labels) {
			for (const auto& item : (*labels))
				output.tag("ResponseType", item.c_str());
		}
		MaterialResponse* resp = new MaterialResponse(this, rid, v);
		output.endTag();
		return resp;
		};

	// the same six labels as the rest of the family
	static std::vector<std::string> lb_tensor = { "11", "22", "33", "12", "23", "13" };
	static std::vector<std::string> lb_ple = { "PLE" };
	static std::vector<std::string> lb_yield = { "F" };
	static std::vector<std::string> lb_error = { "Error" };
	static std::vector<std::string> lb_time = { "dTime", "dTimeCommit", "dTimeInitial" };

	if (argc > 0) {
		// 2000 - model quantities
		if (strcmp(argv[0], "equivalentPlasticStrain") == 0 ||
			strcmp(argv[0], "EquivalentPlasticStrain") == 0 ||
			strcmp(argv[0], "PLE") == 0)
			return make_resp(2000, getEquivalentPlasticStrain(), &lb_ple);
		if (strcmp(argv[0], "backStress1") == 0 || strcmp(argv[0], "BackStress1") == 0)
			return make_resp(2001, getBackStress1(), &lb_tensor);
		if (strcmp(argv[0], "backStress2") == 0 || strcmp(argv[0], "BackStress2") == 0)
			return make_resp(2002, getBackStress2(), &lb_tensor);
		if (strcmp(argv[0], "plasticStrain") == 0 || strcmp(argv[0], "PlasticStrain") == 0)
			return make_resp(2003, getPlasticStrain(), &lb_tensor);
		if (strcmp(argv[0], "yieldFunction") == 0 || strcmp(argv[0], "YieldFunction") == 0)
			return make_resp(2004, getYieldFunction(), &lb_yield);
		// 3000 - implex
		if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0)
			return make_resp(3000, getImplexError(), &lb_error);
		// 3003 - THE STRESS THIS STEP DELIVERED, which under IMPL-EX is not the
		// one 'stress' reports afterwards: commitState installs the implicit
		// solution over it. Same id and same name as the rest of the family -
		// note that ASDSteel1D publishes the same quantity under the same name
		// but at 1009, which is why it is missing from the family-wide test
		if (strcmp(argv[0], "implexStress") == 0 || strcmp(argv[0], "ImplexStress") == 0)
			return make_resp(3003, getImplexStress(), &lb_tensor);
		// 4000 - internal time. Both spellings, because the family does not
		// agree with itself on this one
		if (strcmp(argv[0], "time") == 0 || strcmp(argv[0], "Time") == 0 ||
			strcmp(argv[0], "timeIncrements") == 0 || strcmp(argv[0], "TimeIncrements") == 0)
			return make_resp(4000, getTimeIncrements(), &lb_time);
	}

	// the base class answers "stress", "strain" and "Tangent", which is what the
	// probes and the recorders read - falling through is not optional
	return NDMaterial::setResponse(argv, argc, output);
}

int ASDSteel3DMaterial::getResponse(int responseID, Information& matInformation)
{
	switch (responseID) {
	case 2000: return matInformation.setVector(getEquivalentPlasticStrain());
	case 2001: return matInformation.setVector(getBackStress1());
	case 2002: return matInformation.setVector(getBackStress2());
	case 2003: return matInformation.setVector(getPlasticStrain());
	case 2004: return matInformation.setVector(getYieldFunction());
	case 3000: return matInformation.setVector(getImplexError());
	case 3003: return matInformation.setVector(getImplexStress());
	case 4000: return matInformation.setVector(getTimeIncrements());
	default: break;
	}
	return NDMaterial::getResponse(responseID, matInformation);
}

const Vector& ASDSteel3DMaterial::getEquivalentPlasticStrain() const
{
	out_ple(0) = p;
	return out_ple;
}

const Vector& ASDSteel3DMaterial::getBackStress1() const
{
	out_a1 = a1;
	return out_a1;
}

const Vector& ASDSteel3DMaterial::getBackStress2() const
{
	out_a2 = a2;
	return out_a2;
}

const Vector& ASDSteel3DMaterial::getPlasticStrain() const
{
	// a strain leaves in the ENGINEERING convention, like getStrain()
	for (int i = 0; i < 3; ++i)
		out_ep(i) = ep(i);
	for (int i = 3; i < 6; ++i)
		out_ep(i) = 2.0 * ep(i);
	return out_ep;
}

const Vector& ASDSteel3DMaterial::getYieldFunction() const
{
	static Vector s(6);
	static Vector xi(6);
	deviatorOf(stress, s);
	for (int i = 0; i < 6; ++i)
		xi(i) = s(i) - a1(i) - a2(i);
	out_yield(0) = std::sqrt(std::max(0.0, contract(xi, xi))) - SQ23 * params.sy;
	return out_yield;
}

const Vector& ASDSteel3DMaterial::getImplexError() const
{
	out_error(0) = implex_error;
	return out_error;
}

const Vector& ASDSteel3DMaterial::getImplexStress() const
{
	out_impstress = stress_implex;
	return out_impstress;
}

const Vector& ASDSteel3DMaterial::getTimeIncrements() const
{
	out_time(0) = dtime_n;
	out_time(1) = dtime_n_commit;
	out_time(2) = dtime_0;
	return out_time;
}

// ==================================================================== //
// 12  IMPL-EX error control                                            //
// ==================================================================== //

double ASDSteel3DMaterial::stressReference(void) const
{
	// the largest stress this material can carry. The Chaboche backstresses
	// saturate at H_i/gamma_i, so the asymptote is sy + H1/gamma1 + H2/gamma2 -
	// which, with the calibration the parser does, is exactly the 'su' the user
	// typed. Recovering it from the parameters instead of storing it keeps the
	// serialization alone and stays right if the calibration changes. With no
	// saturation (gamma = 0) the hardening is unbounded and the yield stress is
	// the only scale there is.
	//
	// IDENTICAL to ASDSteel1DMaterial::stressReference(), and that is a
	// requirement, not a coincidence: the numerator of the metric may differ
	// between models (an absolute value in 1D, the largest Voigt component here)
	// but the denominator must not, or the same tolerance means different things
	// in 1D and in 3D and a stepper compares numbers that are not comparable
	double ref = params.sy;
	if (params.gamma1 > 0.0)
		ref += params.H1 / params.gamma1;
	if (params.gamma2 > 0.0)
		ref += params.H2 / params.gamma2;
	return ref > 0.0 ? ref : 1.0;
}

double ASDSteel3DMaterial::implexStressGap(const Vector& delivered, const Vector& stress_implicit) const
{
	double gap = 0.0;
	for (int i = 0; i < 6; ++i)
		gap = std::max(gap, std::abs(delivered(i) - stress_implicit(i)));
	return gap / stressReference();
}

void ASDSteel3DMaterial::saveTrialState(TrialState& x) const
{
	x.stress = stress;
	x.ep = ep;
	x.a1 = a1;
	x.a2 = a2;
	x.n_flow = n_flow;
	x.C = C;
	x.p = p;
	x.dgamma = dgamma;
	x.plastic = plastic;
	x.failed = failed;
}

void ASDSteel3DMaterial::restoreTrialState(const TrialState& x)
{
	stress = x.stress;
	ep = x.ep;
	a1 = x.a1;
	a2 = x.a2;
	n_flow = x.n_flow;
	C = x.C;
	p = x.p;
	dgamma = x.dgamma;
	plastic = x.plastic;
	failed = x.failed;
}

double ASDSteel3DMaterial::computeImplexErrorMetric(void)
{
	// no extrapolation, no error. This is what makes it safe for the aggregation
	// to call it on everything it has
	if (!params.implex)
		return 0.0;

	// the current state holds the EXPLICIT answer: the stress this step
	// delivered to the element, and the state the recorders must keep seeing.
	// stress_implex already holds the delivered stress - setTrialStrain records
	// it - and integrate() does not write that member, so it survives the
	// throw-away solve without being part of the record
	static TrialState delivered;
	saveTrialState(delivered);

	// the implicit answer, at the same trial strain and from the same starting
	// point the commit will use
	if (integrate() < 0) {
		restoreTrialState(delivered);
		// a solve that did not converge cannot say anything about its own error,
		// and no metric is not a small metric: whoever reads this must not
		// accept the step
		return std::numeric_limits<double>::quiet_NaN();
	}
	double err = implexStressGap(stress_implex, stress);

	// undo. Measuring is not allowed to move the state: the step may still be
	// rejected, and the recorders run right after
	restoreTrialState(delivered);
	implex_error = err;
	return err;
}

double ASDSteel3DMaterial::implexTimeRatio(void) const
{
	return dtime_0 > 0.0 ? dtime_n / dtime_0 : 1.0;
}
