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

// Massimo Petracca - ASDEA Software, Italy
//
// Concrete Damaged Plasticity with the damage in the prediction. See the header
// for the formulation, the conventions and what the model does not do.
//
// LAYOUT OF THIS FILE, in the order the model is built up:
//
//   1  Voigt utilities and the coverage counters
//   2  the parser, OPS_ASDPlasticDamageConcrete3DMaterial
//   3  ASDCDPHardeningCurve  - a 1D backbone as a CDP hardening curve
//   4  life-cycle
//   5  the surface and the potential
//   6  the reductions and the damaged elastic law
//   7  the return mapping
//   8  IMPL-EX
//   9  state handling
//  10  the tangent
//  11  output
//  12  serialization

#include <ASDPlasticDamageConcrete3DMaterial.h>
#include <ASDSpectralSplit.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <Element.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <ID.h>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ==================================================================== //
//  1  Voigt utilities and the coverage counters                        //
// ==================================================================== //

// THE COUNTERS ARE OFF BY DEFAULT AND COST NOTHING. They exist because a port
// validated to round-off against its reference on every case is only as strong
// as the coverage behind it: an injected error that changes no output can mean
// either "the port is right there" or "no case ever reaches that line", and only
// a counter tells the two apart. Build with -DASDCDP3D_COUNTERS and call
// ASDPlasticDamageConcrete3DReportCounters() to read them.
#ifdef ASDCDP3D_COUNTERS
namespace {
	struct ASDCDP3DCounters {
		long long integrate = 0;
		long long elastic = 0;
		long long plastic = 0;
		long long iterations = 0;
		long long residuals = 0;
		long long newton_ok = 0;          // den > 0, a Newton direction exists
		long long den_nonpositive = 0;    // den <= 0 or not finite
		long long at_apex = 0;            // ... and the state IS the apex
		long long no_flow = 0;            // m is identically zero: nothing to bracket
		long long damage_return = 0;      // ... so the return ran on kappa_t instead
		long long damage_return_ok = 0;
		long long cycle_bracketed = 0;   // the corrector was going round
		long long bracket_below = 0;      // the corrector overshot: root in (0, lam0)
		long long bracket_called = 0;
		long long bracket_found = 0;
		long long bracket_scan_hit = 0;   // the geometric scan landed on the root
		long long bracket_exhausted = 0;  // 30 bisections were not enough
		long long backtrack = 0;          // a halving was needed
		long long backtrack_exhausted = 0;
		long long gate_by_tol = 0;        // trial accepted because it converged
		long long gate_by_descent = 0;    // ... because |F| went down
		long long gate_by_watchdog = 0;   // ... only by the entry-residual bound
		long long lam_out_of_range = 0;
		long long negative_corrector = 0; // the accumulated multiplier was clamped
		long long cap_binds_kt = 0;
		long long cap_binds_kc = 0;
		long long cap_from_strain = 0;
		long long cap_zero = 0;
		long long stagnated = 0;
		long long max_iter = 0;
		long long best_restored = 0;
		long long degenerate_smax = 0;    // s_max repeated: the averaged projector
		long long macaulay_kink = 0;      // s_max inside the band of the kink
		long long deviator_lost = 0;      // dq/ds returned zero
		long long mises_overflow = 0;
		long long q_tiny_floor = 0;       // the representability floor bit
		long long omega_guard = 0;
		long long extrapolate = 0;
		long long time_factor_clamped = 0;
		long long peek = 0;
		long long eigen_error = 0;
	};
	ASDCDP3DCounters& theCounters() {
		static ASDCDP3DCounters c;
		return c;
	}
}
#define ASDCDP3D_COUNT(X) (++theCounters().X)
void ASDPlasticDamageConcrete3DReportCounters()
{
	const ASDCDP3DCounters& c = theCounters();
#define ASDCDP3D_DUMP(X) opserr << "  " #X " = " << static_cast<double>(c.X) << "\n";
	opserr << "ASDPlasticDamageConcrete3D coverage counters:\n";
	ASDCDP3D_DUMP(integrate) ASDCDP3D_DUMP(elastic) ASDCDP3D_DUMP(plastic)
	ASDCDP3D_DUMP(iterations) ASDCDP3D_DUMP(residuals)
	ASDCDP3D_DUMP(newton_ok) ASDCDP3D_DUMP(den_nonpositive) ASDCDP3D_DUMP(at_apex)
	ASDCDP3D_DUMP(no_flow) ASDCDP3D_DUMP(damage_return) ASDCDP3D_DUMP(damage_return_ok)
	ASDCDP3D_DUMP(cycle_bracketed)
	ASDCDP3D_DUMP(bracket_below)
	ASDCDP3D_DUMP(bracket_called) ASDCDP3D_DUMP(bracket_found)
	ASDCDP3D_DUMP(bracket_scan_hit) ASDCDP3D_DUMP(bracket_exhausted)
	ASDCDP3D_DUMP(backtrack) ASDCDP3D_DUMP(backtrack_exhausted)
	ASDCDP3D_DUMP(gate_by_tol) ASDCDP3D_DUMP(gate_by_descent) ASDCDP3D_DUMP(gate_by_watchdog)
	ASDCDP3D_DUMP(lam_out_of_range) ASDCDP3D_DUMP(negative_corrector)
	ASDCDP3D_DUMP(cap_binds_kt) ASDCDP3D_DUMP(cap_binds_kc)
	ASDCDP3D_DUMP(cap_from_strain) ASDCDP3D_DUMP(cap_zero)
	ASDCDP3D_DUMP(stagnated) ASDCDP3D_DUMP(max_iter) ASDCDP3D_DUMP(best_restored)
	ASDCDP3D_DUMP(degenerate_smax) ASDCDP3D_DUMP(macaulay_kink)
	ASDCDP3D_DUMP(deviator_lost) ASDCDP3D_DUMP(mises_overflow)
	ASDCDP3D_DUMP(q_tiny_floor) ASDCDP3D_DUMP(omega_guard)
	ASDCDP3D_DUMP(extrapolate) ASDCDP3D_DUMP(time_factor_clamped)
	ASDCDP3D_DUMP(peek) ASDCDP3D_DUMP(eigen_error)
#undef ASDCDP3D_DUMP
}
#else
#define ASDCDP3D_COUNT(X) ((void)0)
#endif

namespace {

	enum ErrorCodes {
		EC_Generic = -1,
		EC_IMPLEX_Error_Control = -10,
		EC_Eigen_Error = -1000
	};

	// EVERY 6-VECTOR IN THIS FILE CARRIES THE TENSOR COMPONENTS OF THE SHEAR.
	// The engineering convention lives in setTrialStrain and in the two
	// tangents, and nowhere else - see the header for why the boundary is drawn
	// there. The one consequence to keep in mind is that a full contraction
	// doubles the last three terms, which is what contract() is for.

	inline double macaulay(double x) { return x > 0.0 ? x : 0.0; }

	inline double traceOf(const Vector& a) { return a(0) + a(1) + a(2); }

	// a : b, the full contraction of two symmetric tensors
	inline double contract(const Vector& a, const Vector& b) {
		return a(0) * b(0) + a(1) * b(1) + a(2) * b(2)
			+ 2.0 * (a(3) * b(3) + a(4) * b(4) + a(5) * b(5));
	}

	inline double maxAbs(const Vector& a) {
		double m = 0.0;
		for (int i = 0; i < a.Size(); ++i)
			m = std::max(m, std::abs(a(i)));
		return m;
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

	// q = sqrt(3/2 S:S), non-negative
	inline double misesOf(const Vector& s) {
		double p = traceOf(s) / 3.0;
		double d0 = s(0) - p;
		double d1 = s(1) - p;
		double d2 = s(2) - p;
		return std::sqrt(1.5 * (d0 * d0 + d1 * d1 + d2 * d2
			+ 2.0 * (s(3) * s(3) + s(4) * s(4) + s(5) * s(5))));
	}

	// sigma = lam*tr(e)*I + 2*mu*e, in closed form. mu2 = 2*mu.
	//
	// This is the whole dividend of the tensor-component convention: the elastic
	// law needs no 6x6 matrix and no engineering factor, in either direction, and
	// it applies unchanged to a strain and to a flow direction.
	inline void elasticStress(double lam, double mu2, const Vector& e, Vector& out) {
		double lam2 = lam + mu2;
		out(0) = lam2 * e(0) + lam * e(1) + lam * e(2);
		out(1) = lam * e(0) + lam2 * e(1) + lam * e(2);
		out(2) = lam * e(0) + lam * e(1) + lam2 * e(2);
		out(3) = mu2 * e(3);
		out(4) = mu2 * e(4);
		out(5) = mu2 * e(5);
	}

	// The 6x6 elastic operator in the SAME convention as the split projectors:
	// it maps a 6-vector of tensor components to one of tensor components, so its
	// shear columns carry the factor two of the contraction. That is what makes
	// the composition W:C a plain matrix product.
	inline void elasticMatrix(double lam, double mu2, Matrix& out) {
		out.Zero();
		double lam2 = lam + mu2;
		out(0, 0) = out(1, 1) = out(2, 2) = lam2;
		out(0, 1) = out(1, 0) = out(0, 2) = out(2, 0) = out(1, 2) = out(2, 1) = lam;
		out(3, 3) = out(4, 4) = out(5, 5) = mu2;
	}

	// Lee and Fenves' r = sum<s_i> / sum|s_i|: HOW TENSILE a state is, in [0, 1].
	// 1 where the principal stresses are all non-negative, 0 where they are all
	// non-positive, the ratio of the two in between.
	//
	// At s = 0 it is 0/0 - a known annoyance of the CDP - and the choice here is
	// 0.5. It is never reached in practice with the flow rule in front of it,
	// there being no plastic flow at zero stress, and a symmetric value is the
	// only one that does not bias a side.
	//
	// The eigenvalues arrive DESCENDING and are summed from the smallest up,
	// which is the order the reference implementation's ascending array gives.
	inline double positiveFraction(const Vector& d, double tol = 1.0e-14) {
		double tot = std::abs(d(2)) + std::abs(d(1)) + std::abs(d(0));
		if (tot <= tol)
			return 0.5;
		double pos = macaulay(d(2)) + macaulay(d(1)) + macaulay(d(0));
		return pos / tot;
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

	bool string_to_list_of_doubles(const std::string& text, char sep, std::vector<double>& out) {
		if (out.size() > 0) out.clear();
		std::size_t start = 0, end = 0;
		double value;
		while (true) {
			end = text.find(sep, start);
			if (end == std::string::npos) {
				if (start < text.size()) {
					if (!string_to_double(text.substr(start), value))
						return false;
					out.push_back(value);
				}
				break;
			}
			std::string subs = text.substr(start, end - start);
			if (subs.size() > 0) {
				if (!string_to_double(subs, value))
					return false;
				out.push_back(value);
			}
			start = end + 1;
		}
		return true;
	}

}

const double ASDPlasticDamageConcrete3DMaterial::ResidualGrowth = 1.0;

// ==================================================================== //
//  2  the parser                                                       //
// ==================================================================== //

void* OPS_ASDPlasticDamageConcrete3DMaterial(void)
{
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDPlasticDamageConcrete3D - Developed by: Massimo Petracca, ASDEA Software Technology\n";
		first_done = true;
	}

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 3) {
		opserr <<
			"nDMaterial ASDPlasticDamageConcrete3D Error: Few arguments (< 3).\n"
			"nDMaterial ASDPlasticDamageConcrete3D $tag $E $nu "
			"-Te $Te -Ts $Ts -Ce $Ce -Cs $Cs "
			"<-rho $rho> <-dilatancy $psiDegrees> <-eccentricity $ecc> "
			"<-fb0fc0 $ratio> <-Kc $Kc> <-damageT $dt> <-damageC $dc> "
			"<-damageCombination faria|leeFenves> "
			"<-stiffnessRecoveryT $wt> <-stiffnessRecoveryC $wc> "
			"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> "
			"<-implexAbort> <-implexAlpha $alpha> "
			"<-autoRegularization $lch_ref> "
			"<-integration $tol $maxIter $maxBacktrack $stagnationTol>\n";
		return nullptr;
	}

	int numData = 1;

	int tag;
	double E;
	double nu = 0.2;
	double rho = 0.0;
	double dilatancy = 30.0;
	double eccentricity = 0.1;
	double fb0_fc0 = 1.16;
	double Kc = 2.0 / 3.0;
	double damage_t = 1.0;
	double damage_c = 0.0;
	ASDPlasticDamageConcrete3DMaterial::DamageCombination damage_combination =
		ASDPlasticDamageConcrete3DMaterial::DC_Faria;
	double stiffness_recovery_t = 0.0;
	double stiffness_recovery_c = 1.0;
	// WHETHER the two weights were given, not only their value, because a
	// weight that has nothing to weigh is refused rather than ignored - and the
	// options may arrive in any order, so the check cannot live in the loop
	bool recovery_t_given = false;
	bool recovery_c_given = false;
	bool implex = false;
	bool implex_control = false;
	bool implex_abort_on_error = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	double implex_alpha = 1.0;
	bool auto_regularization = false;
	double lch_ref = 1.0;
	double tol = 1.0e-10;
	int max_iter = 100;
	int max_backtrack = 12;
	double stagnation_tol = 1.0e-6;
	std::vector<double> Te, Ts, Ce, Cs;

	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: invalid 'tag'.\n";
		return nullptr;
	}
	if (OPS_GetDouble(&numData, &E) != 0) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: invalid 'E'.\n";
		return nullptr;
	}
	if (E <= 0.0) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: invalid value for 'E' (" << E << "). It should be strictly positive.\n";
		return nullptr;
	}
	if (OPS_GetDouble(&numData, &nu) != 0) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: invalid 'nu'.\n";
		return nullptr;
	}

	auto lam_optional_int = [&numData](const char* variable, int& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetInt(&numData, &value) < 0) {
				opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_list = [&numData](const char* variable, std::vector<double>& value) -> bool {
		// first try an expanded list, {*}$the_list in Tcl or *the_list in python
		value.clear();
		while (OPS_GetNumRemainingInputArgs() > 0) {
			double item;
			auto old_num_rem = OPS_GetNumRemainingInputArgs();
			if (OPS_GetDoubleInput(&numData, &item) < 0) {
				auto new_num_rem = OPS_GetNumRemainingInputArgs();
				if (new_num_rem < old_num_rem)
					OPS_ResetCurrentInputArg(-1);
				break;
			}
			value.push_back(item);
		}
		// then a Tcl list (it is a string after all)
		if (value.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
			std::string list_string = OPS_GetString();
			if (!string_to_list_of_doubles(list_string, ' ', value)) {
				opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: cannot parse the '" << variable << "' list.\n";
				return false;
			}
		}
		return true;
	};

	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-rho") == 0) {
			if (!lam_optional_double("rho", rho))
				return nullptr;
		}
		else if (strcmp(value, "-dilatancy") == 0) {
			if (!lam_optional_double("dilatancy", dilatancy))
				return nullptr;
		}
		else if (strcmp(value, "-eccentricity") == 0) {
			if (!lam_optional_double("eccentricity", eccentricity))
				return nullptr;
		}
		else if (strcmp(value, "-fb0fc0") == 0) {
			if (!lam_optional_double("fb0fc0", fb0_fc0))
				return nullptr;
		}
		else if (strcmp(value, "-Kc") == 0) {
			if (!lam_optional_double("Kc", Kc))
				return nullptr;
		}
		else if (strcmp(value, "-damageT") == 0) {
			if (!lam_optional_double("damageT", damage_t))
				return nullptr;
		}
		else if (strcmp(value, "-damageC") == 0) {
			if (!lam_optional_double("damageC", damage_c))
				return nullptr;
		}
		else if (strcmp(value, "-damageCombination") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: '-damageCombination' given without the next 1 argument (faria|leeFenves).\n";
				return nullptr;
			}
			const char* how = OPS_GetString();
			if (strcmp(how, "faria") == 0 || strcmp(how, "Faria") == 0)
				damage_combination = ASDPlasticDamageConcrete3DMaterial::DC_Faria;
			else if (strcmp(how, "leeFenves") == 0 || strcmp(how, "LeeFenves") == 0 ||
				strcmp(how, "leefenves") == 0)
				damage_combination = ASDPlasticDamageConcrete3DMaterial::DC_LeeFenves;
			else {
				opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: unknown '-damageCombination' value '" <<
					how << "'. It is 'faria' - the two reductions applied to the two "
					"spectral parts separately, where the recovery of stiffness on a "
					"reversal is total and automatic - or 'leeFenves' - the two "
					"combined into the single scalar of the classical CDP, where the "
					"recovery is governed by '-stiffnessRecoveryT'/'-stiffnessRecoveryC'.\n";
				return nullptr;
			}
		}
		else if (strcmp(value, "-stiffnessRecoveryT") == 0) {
			if (!lam_optional_double("stiffnessRecoveryT", stiffness_recovery_t))
				return nullptr;
			recovery_t_given = true;
		}
		else if (strcmp(value, "-stiffnessRecoveryC") == 0) {
			if (!lam_optional_double("stiffnessRecoveryC", stiffness_recovery_c))
				return nullptr;
			recovery_c_given = true;
		}
		else if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
				return nullptr;
			}
			if (!lam_optional_double("implexErrorTolerance", implex_error_tolerance))
				return nullptr;
			if (!lam_optional_double("implexTimeReductionLimit", implex_time_redution_limit))
				return nullptr;
		}
		else if (strcmp(value, "-implexAbort") == 0) {
			// LEGACY: let the material fail the step by itself. See the note on
			// implex_abort_on_error in the header - the rejection belongs to
			// CTestImplexWrapper
			implex_abort_on_error = true;
		}
		else if (strcmp(value, "-implexAlpha") == 0) {
			if (!lam_optional_double("alpha", implex_alpha))
				return nullptr;
		}
		else if (strcmp(value, "-autoRegularization") == 0) {
			auto_regularization = true;
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: '-autoRegularization' given without the next 1 argument $lch_ref.\n";
				return nullptr;
			}
			if (!lam_optional_double("lch_ref", lch_ref))
				return nullptr;
		}
		else if (strcmp(value, "-integration") == 0) {
			if (OPS_GetNumRemainingInputArgs() < 4) {
				opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: '-integration' given without the next 4 arguments $tol $maxIter $maxBacktrack $stagnationTol.\n";
				return nullptr;
			}
			if (!lam_optional_double("tol", tol))
				return nullptr;
			if (!lam_optional_int("maxIter", max_iter))
				return nullptr;
			if (!lam_optional_int("maxBacktrack", max_backtrack))
				return nullptr;
			if (!lam_optional_double("stagnationTol", stagnation_tol))
				return nullptr;
		}
		else if (strcmp(value, "-Te") == 0) {
			if (!lam_optional_list("Te", Te))
				return nullptr;
		}
		else if (strcmp(value, "-Ts") == 0) {
			if (!lam_optional_list("Ts", Ts))
				return nullptr;
		}
		else if (strcmp(value, "-Ce") == 0) {
			if (!lam_optional_list("Ce", Ce))
				return nullptr;
		}
		else if (strcmp(value, "-Cs") == 0) {
			if (!lam_optional_list("Cs", Cs))
				return nullptr;
		}
		else if (strcmp(value, "-Td") == 0 || strcmp(value, "-Cd") == 0) {
			// REFUSED, NOT IGNORED, and that is the point. The tensile reduction
			// of this model is DERIVED from the nominal backbone,
			// omega_t = q/(df*E*kappa + q), and compression at the default dial
			// is pure plasticity: a damage array has no place to go here. It
			// would not even be inert - ASDHardeningLaw::adjust() reshapes the
			// backbone from it - so accepting and dropping it would change the
			// curve the user thinks they gave. A silently dropped input is worse
			// than a rejected one
			opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: '" << value <<
				"' is not an input of this model: the tensile reduction is derived "
				"from the NOMINAL backbone (omega_t = q/(df*E*kappa+q)) and the "
				"compressive side is pure plasticity at the default dial. Pass the "
				"nominal curves only, and move the split between damage and "
				"plasticity onto '-damageT'/'-damageC'.\n";
			return nullptr;
		}
		else if (strcmp(value, "-fc") == 0 || strcmp(value, "-ft") == 0) {
			// NO GENERATED PRESET HERE, on purpose. The generator of
			// ASDConcrete3D produces the damage arrays along with the two
			// backbones, and this model refuses those (see above). Reusing only
			// its Te/Ts/Ce/Cs would hand over a curve whose shape was calibrated
			// against a damage split this model derives differently, which is
			// exactly the kind of quiet mismatch that is worth one explicit
			// input instead
			opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: '" << value <<
				"' is not an input of this model: the generated preset of "
				"ASDConcrete3D also generates the damage arrays, which this model "
				"derives instead. Pass the two nominal backbones explicitly with "
				"-Te/-Ts and -Ce/-Cs.\n";
			return nullptr;
		}
	}

	// --- bounds. EVERY ONE OF THESE SITS IN A DENOMINATOR of a closed form
	// below, so the failure without them is a division by zero three frames deep
	// instead of a rejected input. The bounds are not taste, they are where the
	// formulas stop existing: fb0/fc0 = 0.5 kills alpha, Kc = 0.5 kills gamma,
	// nu = 0.5 kills the bulk modulus and psi = 90 deg kills tan(psi)
	// ... AND FOR nu THE BOUND IS TIGHTER THAN WHERE THE FORMULA DIES, on purpose.
	// Thermodynamics allows (-1, 0.5) and the elastic tensor exists on all of it,
	// but this is a CONCRETE model: nu < 0 is auxetic, which concrete is not, and
	// the negative end was not merely unused - it was actively misleading. The
	// 55-case validation suite carried a nu = -0.5 case, and that one case
	// produced 43 of the 380 failed steps, 43 of the 135 exhausted backtracks and
	// the WORST disagreement with the Python bench in the whole suite, 95.6 MPa.
	// Every summary of the model's robustness was dominated by an input nobody
	// will ever write.
	//
	// The upper end is 0.499 and not 0.5 for the same practical reason: the pole
	// is at 0.5, but lam = nu*mu2/(1-2nu) is already 249*mu2 at 0.499 and 2.5e6
	// times mu2 at 0.4999999, so 'representable' stops meaning 'usable' well
	// before the formula stops existing. Anything closer to the pole is a
	// conditioning problem the caller cannot see, so it is refused where it can.
	if (!(nu >= 0.0 && nu <= 0.499)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'nu' must be in [0, 0.499], got " << nu <<
			": below 0 is auxetic, which concrete is not, and above 0.499 the bulk"
			" modulus is so close to its pole at 0.5 that lam = nu*mu2/(1-2nu)"
			" swamps the shear modulus.\n";
		return nullptr;
	}
	if (!(fb0_fc0 > 0.5)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'fb0fc0' must be > 0.5, got " << fb0_fc0 <<
			": alpha = (r-1)/(2r-1) is what it enters, and 0.5 is its pole.\n";
		return nullptr;
	}
	if (!(Kc > 0.5 && Kc <= 1.0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'Kc' must be in (0.5, 1], got " << Kc <<
			": gamma = 3(1-Kc)/(2Kc-1) is what it enters, and 0.5 is its pole.\n";
		return nullptr;
	}
	if (!(dilatancy >= 0.0 && dilatancy < 90.0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'dilatancy' must be in [0, 90) degrees, got " << dilatancy << ".\n";
		return nullptr;
	}
	if (!(eccentricity >= 0.0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'eccentricity' must be >= 0, got " << eccentricity << ".\n";
		return nullptr;
	}
	if (!(damage_t >= 0.0 && damage_t <= 1.0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'damageT' must be in [0, 1], got " << damage_t << ".\n";
		return nullptr;
	}
	if (!(damage_c >= 0.0 && damage_c <= 1.0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'damageC' must be in [0, 1], got " << damage_c << ".\n";
		return nullptr;
	}
	// THE TWO WEIGHTS ARE REFUSED WHERE THEY WOULD DO NOTHING, not accepted and
	// dropped, which is the same position this parser takes on '-Td'/'-Cd' and
	// on '-fc'. In the Faria split the recovery of stiffness is performed by the
	// spectral split itself and is total: a cracked side loses only its positive
	// part, so there is no fraction left for a weight to choose. Accepting the
	// input would tell the user they had calibrated something they had not
	if ((recovery_t_given || recovery_c_given) &&
		damage_combination != ASDPlasticDamageConcrete3DMaterial::DC_LeeFenves) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: '-stiffnessRecoveryT'"
			"/'-stiffnessRecoveryC' exist only with '-damageCombination leeFenves'. In"
			" the Faria split the two reductions are applied to the two spectral parts"
			" separately, so the stiffness of a side is recovered in full and"
			" automatically as soon as the stress changes sign, and there is no"
			" fraction left for a weight to govern.\n";
		return nullptr;
	}
	if (!(stiffness_recovery_t >= 0.0 && stiffness_recovery_t <= 1.0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'stiffnessRecoveryT' must be in [0, 1], got " <<
			stiffness_recovery_t << ": it is the fraction of the TENSILE stiffness"
			" recovered from the compressive damage when the state turns tensile, and"
			" outside [0, 1] the scalar reduction leaves (0, 1].\n";
		return nullptr;
	}
	if (!(stiffness_recovery_c >= 0.0 && stiffness_recovery_c <= 1.0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'stiffnessRecoveryC' must be in [0, 1], got " <<
			stiffness_recovery_c << ": it is the fraction of the COMPRESSIVE stiffness"
			" recovered from the tensile damage when the crack closes, and outside"
			" [0, 1] the scalar reduction leaves (0, 1].\n";
		return nullptr;
	}
	if (!(tol > 0.0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'tol' must be > 0, got " << tol << ".\n";
		return nullptr;
	}
	if (!(max_iter >= 1)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'maxIter' must be >= 1, got " << max_iter << ".\n";
		return nullptr;
	}
	if (!(max_backtrack >= 0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'maxBacktrack' must be >= 0, got " << max_backtrack << ".\n";
		return nullptr;
	}
	if (!(stagnation_tol >= 0.0)) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'stagnationTol' must be >= 0, got " << stagnation_tol << ".\n";
		return nullptr;
	}

	// --- the two backbones
	if (Te.size() < 1) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'Te' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Ts.size() != Te.size()) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Ts' (size = " <<
			static_cast<int>(Ts.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Ce.size() < 1) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'Ce' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Cs.size() != Ce.size()) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cs' (size = " <<
			static_cast<int>(Cs.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	// the damage arrays are ZERO, by construction and not by omission: this
	// model's reduction is derived from the curve
	std::vector<double> Td(Te.size(), 0.0);
	std::vector<double> Cd(Ce.size(), 0.0);

	ASDPlasticDamageConcrete3DMaterial::HardeningLaw HT(
		tag, ASDHardeningLawType::Tension, E, Te, Ts, Td);
	if (!HT.isValid()) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: Tensile hardening law is not valid.\n";
		return nullptr;
	}
	ASDPlasticDamageConcrete3DMaterial::HardeningLaw HC(
		tag, ASDHardeningLawType::Compression, E, Ce, Cs, Cd);
	if (!HC.isValid()) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: Compressive hardening law is not valid.\n";
		return nullptr;
	}

	// A DRY RUN OF THE CONVERSION, HERE, where the user can still act on it. The
	// real one is lazy and has to be: it happens after regularization, which
	// needs the parent element's characteristic length and is therefore only
	// possible at the first setTrialStrain. A table that cannot become a CDP
	// hardening curve would then be reported from inside an analysis, which is
	// the wrong place to learn it
	{
		ASDCDPHardeningCurve probe;
		const char* why = 0;
		if (!probe.build(HT, E, 0.0, &why)) {
			opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: the tensile backbone cannot be used as a hardening curve: " << why << "\n";
			return nullptr;
		}
		if (!probe.build(HC, E, 0.0, &why)) {
			opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: the compressive backbone cannot be used as a hardening curve: " << why << "\n";
			return nullptr;
		}
	}

	NDMaterial* instance = new ASDPlasticDamageConcrete3DMaterial(
		tag,
		E, nu, rho,
		dilatancy, eccentricity, fb0_fc0, Kc,
		damage_t, damage_c,
		damage_combination, stiffness_recovery_t, stiffness_recovery_c,
		implex, implex_control, implex_abort_on_error,
		implex_error_tolerance, implex_time_redution_limit, implex_alpha,
		auto_regularization, lch_ref,
		HT, HC,
		tol, max_iter, max_backtrack, stagnation_tol);
	if (instance == nullptr) {
		opserr << "nDMaterial ASDPlasticDamageConcrete3D Error: failed to allocate a new material.\n";
		return nullptr;
	}
	return instance;
}

// ==================================================================== //
//  3  ASDCDPHardeningCurve                                             //
// ==================================================================== //

bool ASDCDPHardeningCurve::build(const ASDHardeningLaw& law, double E,
	double damage_cap, const char** why)
{
	m_kappa.clear();
	m_q.clear();
	m_y.clear();
	m_qmax = 0.0;
	m_dropped = 0;
	if (why) *why = "";

	if (!law.isValid() || law.points().size() < 2) {
		if (why) *why = "the hardening law is not valid";
		return false;
	}
	if (!(E > 0.0)) {
		if (why) *why = "E must be strictly positive";
		return false;
	}

	const std::vector<ASDHardeningLawPoint>& pts = law.points();
	double xtol = law.strainTolerance();
	// THE ORIGIN IS DROPPED: the curve starts at kappa = 0 carrying the ELASTIC
	// LIMIT, which is where a yield surface starts
	for (std::size_t i = 1; i < pts.size(); ++i) {
		const ASDHardeningLawPoint& p = pts[i];
		// CAP THE DAMAGE, because it is what the effective strength is divided
		// by. At 0 - this model's case - the conversion is the identity in
		// stress and the curve IS the nominal backbone. Uncapped, a Model Code
		// tensile array ends at d = 0.999993 and asks for q = 434.85 MPa on a
		// 3 MPa material, with two adjacent points a hardening modulus of 4e20
		// apart
		double d = std::min(p.d, damage_cap);
		double q = (d < 1.0) ? p.y / (1.0 - d) : p.q;
		double k = p.x - q / E;
		if (m_kappa.size() > 0 && k <= m_kappa.back() + xtol) {
			// keep the FIRST value at a repeated kappa: the curve may not jump,
			// and the later point is the one asking it to
			++m_dropped;
			continue;
		}
		m_kappa.push_back(k);
		m_q.push_back(q);
		m_y.push_back(p.y);
	}
	if (m_kappa.size() == 0) {
		if (why) *why = "no usable hardening point: no point advances the hardening measure kappa = x - q/E";
		return false;
	}
	// A NON-POSITIVE EFFECTIVE STRENGTH IS NOT A CURVE THIS SURFACE CAN USE, and
	// saying so here is worth more than any guard downstream: beta divides by qt
	// and dF/dkt by qt*qt, so a zero would come out as a division by zero from
	// inside the return mapping. ASDHardeningLaw already clamps y to its stress
	// tolerance, so this fires only on a law built by hand or mutated after the
	// fact
	for (std::size_t i = 0; i < m_q.size(); ++i) {
		if (!(m_q[i] > 0.0) || !std::isfinite(m_q[i])) {
			if (why) *why = "the effective strength must be > 0 and finite at every point: the yield surface divides by it";
			return false;
		}
	}
	// the elastic limit sits at kappa = 0 by construction (d = 0 at the first
	// corner is enforced by adjust()), but say so rather than assume it
	if (m_kappa[0] > xtol) {
		m_kappa.insert(m_kappa.begin(), 0.0);
		m_q.insert(m_q.begin(), m_q[0]);
		m_y.insert(m_y.begin(), m_y[0]);
	}
	for (std::size_t i = 0; i < m_q.size(); ++i)
		m_qmax = std::max(m_qmax, m_q[i]);
	return true;
}

void ASDCDPHardeningCurve::evaluate(double kappa, double& q, double& dq) const
{
	std::size_t n = m_kappa.size();
	if (n == 0) {
		q = 0.0;
		dq = 0.0;
		return;
	}
	if (kappa <= m_kappa[0]) {
		// AT the first point the slope must be the FORWARD one, not zero - see
		// the header
		q = m_q[0];
		dq = (n > 1) ? (m_q[1] - m_q[0]) / (m_kappa[1] - m_kappa[0]) : 0.0;
		return;
	}
	if (kappa >= m_kappa[n - 1]) {
		// held CONSTANT and not extrapolated along the last tangent: a negative
		// one would drive the effective strength through zero and take the yield
		// surface with it
		q = m_q[n - 1];
		dq = 0.0;
		return;
	}
	// the largest i with kappa[i] <= kappa, which is what searchsorted(...,
	// 'right') - 1 gives
	std::size_t i = 0;
	for (std::size_t j = 1; j < n; ++j) {
		if (m_kappa[j] <= kappa)
			i = j;
		else
			break;
	}
	double h = (m_q[i + 1] - m_q[i]) / (m_kappa[i + 1] - m_kappa[i]);
	q = m_q[i] + h * (kappa - m_kappa[i]);
	dq = h;
}

double ASDCDPHardeningCurve::nominal(double kappa) const
{
	std::size_t n = m_kappa.size();
	if (n == 0)
		return 0.0;
	if (kappa <= m_kappa[0])
		return m_y[0];
	if (kappa >= m_kappa[n - 1])
		return m_y[n - 1];
	std::size_t i = 0;
	for (std::size_t j = 1; j < n; ++j) {
		if (m_kappa[j] <= kappa)
			i = j;
		else
			break;
	}
	double t = (kappa - m_kappa[i]) / (m_kappa[i + 1] - m_kappa[i]);
	return m_y[i] + t * (m_y[i + 1] - m_y[i]);
}

double ASDCDPHardeningCurve::omegaOfCurve(double kappa) const
{
	double q, dq;
	evaluate(kappa, q, dq);
	return (q > 0.0) ? nominal(kappa) / q : 1.0;
}

// ==================================================================== //
//  4  life-cycle                                                       //
// ==================================================================== //

ASDPlasticDamageConcrete3DMaterial::ASDPlasticDamageConcrete3DMaterial(
	int _tag,
	double _E,
	double _nu,
	double _rho,
	double _dilatancy_deg,
	double _eccentricity,
	double _fb0_fc0,
	double _Kc,
	double _damage_t,
	double _damage_c,
	DamageCombination _damage_combination,
	double _stiffness_recovery_t,
	double _stiffness_recovery_c,
	bool _implex,
	bool _implex_control,
	bool _implex_abort_on_error,
	double _implex_error_tolerance,
	double _implex_time_reduction_limit,
	double _implex_alpha,
	bool _auto_regularize,
	double _lch_ref,
	const HardeningLaw& _ht,
	const HardeningLaw& _hc,
	double _tol,
	int _max_iter,
	int _max_backtrack,
	double _stagnation_tol)
	: NDMaterial(_tag, ND_TAG_ASDPlasticDamageConcrete3DMaterial)
	, E(_E)
	, nu(_nu)
	, rho(_rho)
	, psi(std::abs(_dilatancy_deg)* M_PI / 180.0)
	, ecc(_eccentricity)
	, fb0_fc0(_fb0_fc0)
	, Kc(_Kc)
	, damage_t(_damage_t)
	, damage_c(_damage_c)
	, damage_combination(_damage_combination)
	, stiffness_recovery_t(_stiffness_recovery_t)
	, stiffness_recovery_c(_stiffness_recovery_c)
	, implex(_implex)
	, implex_control(_implex_control)
	, implex_abort_on_error(_implex_abort_on_error)
	, implex_error_tolerance(_implex_error_tolerance)
	, implex_time_redution_limit(_implex_time_reduction_limit)
	, implex_alpha(_implex_alpha)
	, tol(_tol)
	, max_iter(_max_iter)
	, max_backtrack(_max_backtrack)
	, stagnation_tol(_stagnation_tol)
	, auto_regularize(_auto_regularize)
	, lch_ref(_lch_ref)
	, ht(_ht)
	, hc(_hc)
{
	// Lubliner's alpha, from the equibiaxial ratio, and the third-invariant
	// shape factor, active only where s_max < 0
	alpha = (fb0_fc0 - 1.0) / (2.0 * fb0_fc0 - 1.0);
	gam = 3.0 * (1.0 - Kc) / (2.0 * Kc - 1.0);
	// the operator is the identity while both reductions are 1, and the frozen
	// split starts EVEN: zero eigenvalues put the whole projector in the shared
	// remainder, so PT = PC = I/2, which is what the reference implementation
	// initializes them to explicitly
	W.Zero();
	for (int i = 0; i < 6; ++i)
		W(i, i) = 1.0;
	V_split.Zero();
	V_split_commit.Zero();
	for (int i = 0; i < 3; ++i) {
		V_split(i, i) = 1.0;
		V_split_commit(i, i) = 1.0;
	}
	C = getInitialTangent();
}

ASDPlasticDamageConcrete3DMaterial::ASDPlasticDamageConcrete3DMaterial()
	: NDMaterial(0, ND_TAG_ASDPlasticDamageConcrete3DMaterial)
{
}

ASDPlasticDamageConcrete3DMaterial::~ASDPlasticDamageConcrete3DMaterial()
{
}

double ASDPlasticDamageConcrete3DMaterial::getRho(void)
{
	return rho;
}

NDMaterial* ASDPlasticDamageConcrete3DMaterial::getCopy(void)
{
	// the default copy-constructor is safe for the members this class uses
	return new ASDPlasticDamageConcrete3DMaterial(*this);
}

NDMaterial* ASDPlasticDamageConcrete3DMaterial::getCopy(const char* code)
{
	if (strcmp(code, "ThreeDimensional") == 0)
		return getCopy();
	return NDMaterial::getCopy(code);
}

const char* ASDPlasticDamageConcrete3DMaterial::getType(void) const
{
	return "ThreeDimensional";
}

int ASDPlasticDamageConcrete3DMaterial::getOrder(void) const
{
	return 6;
}

bool ASDPlasticDamageConcrete3DMaterial::prepare(void)
{
	if (curves_ready)
		return setup_ok;

	// REGULARIZE ONCE, and only on the build path: on the receive path the laws
	// arrive already regularized and regularization_done travels with them
	if (!regularization_done) {
		if (ops_TheActiveElement)
			lch = ops_TheActiveElement->getCharacteristicLength();
		regularization_done = true;
		if (auto_regularize) {
			ht.regularize(lch, lch_ref);
			hc.regularize(lch, lch_ref);
		}
	}

	curves_ready = true;
	const char* why_t = 0;
	const char* why_c = 0;
	// THE DAMAGE CAP IS ZERO, which is what makes these the NOMINAL backbones
	// (abscissa x - y/E, ordinate y) and the input damage arrays irrelevant: the
	// reduction of this model is derived from the curve itself
	bool ok_t = hard_t.build(ht, E, 0.0, &why_t);
	bool ok_c = hard_c.build(hc, E, 0.0, &why_c);
	setup_ok = ok_t && ok_c;
	if (!setup_ok) {
		if (!setup_reported) {
			setup_reported = true;
			opserr << "ASDPlasticDamageConcrete3DMaterial (tag = " << getTag() <<
				") Error: the backbones cannot be converted into hardening curves.\n";
			if (!ok_t)
				opserr << "   tension: " << why_t << "\n";
			if (!ok_c)
				opserr << "   compression: " << why_c << "\n";
		}
		return false;
	}
	if (!setup_reported && (hard_t.dropped() > 0 || hard_c.dropped() > 0)) {
		// NOT SILENT. A point that does not advance kappa is not interpolable -
		// q jumping at constant kappa is an infinite hardening modulus - and
		// smoothing over it would hide the one thing the user needs to be told
		// about their input
		setup_reported = true;
		opserr << "ASDPlasticDamageConcrete3DMaterial (tag = " << getTag() <<
			") Warning: hardening points that do not advance kappa = x - q/E were "
			"dropped (tension: " << hard_t.dropped() << ", compression: " <<
			hard_c.dropped() << "). They are not interpolable, and the value kept "
			"at a repeated kappa is the FIRST one.\n";
	}
	// 24 orders of magnitude below the strongest thing the material carries, and
	// never below the point where qt*qt stops being a number
	double big = std::max(std::max(hard_t.maxEffectiveStress(),
		hard_c.maxEffectiveStress()), 1.0);
	q_tiny = std::max(1.0e-24 * big, 1.0e-150);
	return true;
}

// ==================================================================== //
//  5  the surface and the potential                                    //
// ==================================================================== //

double ASDPlasticDamageConcrete3DMaterial::surfaceBeta(double qt, double qc) const
{
	double qtf = std::max(qt, q_tiny);
#ifdef ASDCDP3D_COUNTERS
	if (qtf != qt) ASDCDP3D_COUNT(q_tiny_floor);
#endif
	return qc / qtf * (1.0 - alpha) - (1.0 + alpha);
}

double ASDPlasticDamageConcrete3DMaterial::yieldFunction(const Vector& s,
	double qt, double qc) const
{
	static Vector d(3);
	if (ASDSpectralSplit::eigenvalues(s, d) < 0) {
		ASDCDP3D_COUNT(eigen_error);
		return std::numeric_limits<double>::quiet_NaN();
	}
	double smax = d(0);
	double beta = surfaceBeta(qt, qc);
	double val = misesOf(s) + alpha * traceOf(s)
		+ beta * macaulay(smax) - gam * macaulay(-smax);
	return val / (1.0 - alpha) - qc;
}

void ASDPlasticDamageConcrete3DMaterial::dFdKappa(double smax,
	double qt, double qc, double dqt, double dqc,
	double& dF_dkt, double& dF_dkc) const
{
	double qtf = std::max(qt, q_tiny);
	double sm = macaulay(smax);
	dF_dkt = sm * (-qc * (1.0 - alpha) * dqt / (qtf * qtf)) / (1.0 - alpha);
	dF_dkc = sm * (dqc * (1.0 - alpha) / qtf) / (1.0 - alpha) - dqc;
}

void ASDPlasticDamageConcrete3DMaterial::smaxProjector(const Vector& s,
	double& smax, Vector& P) const
{
	static Vector d(3);
	static Matrix V(3, 3);
	P.Zero();
	if (ASDSpectralSplit::spectral(s, d, V) < 0) {
		ASDCDP3D_COUNT(eigen_error);
		smax = 0.0;
		return;
	}
	smax = d(0);
	double amax = std::max(std::max(std::abs(d(0)), std::abs(d(1))), std::abs(d(2)));
	double scale = std::max(amax, 1.0);
	double etol = 1.0e-10 * scale;
	int count = 0;
	for (int i = 0; i < 3; ++i) {
		if (d(i) >= smax - etol) {
			++count;
			// n (x) n of that direction, in tensor components
			double nx = V(0, i), ny = V(1, i), nz = V(2, i);
			P(0) += nx * nx;
			P(1) += ny * ny;
			P(2) += nz * nz;
			P(3) += nx * ny;
			P(4) += ny * nz;
			P(5) += nx * nz;
		}
	}
#ifdef ASDCDP3D_COUNTERS
	if (count > 1) ASDCDP3D_COUNT(degenerate_smax);
#endif
	if (count > 1) {
		double f = 1.0 / static_cast<double>(count);
		for (int i = 0; i < 6; ++i)
			P(i) *= f;
	}
}

void ASDPlasticDamageConcrete3DMaterial::dqds(const Vector& s, Vector& out) const
{
	static Vector S(6);
	deviatorOf(s, S);
	double scale = maxAbs(S);
	if (!(scale > 0.0) || !std::isfinite(scale)) {
		ASDCDP3D_COUNT(deviator_lost);
		out.Zero();
		return;
	}
	// SCALE-INVARIANT BY CONSTRUCTION - it is sqrt(3/2) times the unit deviator -
	// so normalizing changes no value and removes the only way the expression
	// can fail
	for (int i = 0; i < 6; ++i)
		S(i) /= scale;
	double f = 1.5 / std::sqrt(1.5 * contract(S, S));
	for (int i = 0; i < 6; ++i)
		out(i) = f * S(i);
}

void ASDPlasticDamageConcrete3DMaterial::dFds(const Vector& s, double qt, double qc,
	Vector& out, double& smax) const
{
	static Vector P(6);
	static Vector dq(6);
	double beta = surfaceBeta(qt, qc);
	smaxProjector(s, smax, P);
	dqds(s, dq);
	// The Macaulay pair has a KINK at s_max = 0, from slope gamma below to slope
	// beta above. Both are positive and, for concrete, close (3.00 and 3.76 at
	// first yield), so the surface is nearly smooth there - but uniaxial AND
	// equibiaxial compression both sit exactly ON it, so the value at zero
	// cannot be left to a strict inequality. The midpoint is the symmetric
	// element of the subdifferential, the same choice smaxProjector makes
	double stol = 1.0e-12 * std::max(qc, 1.0);
	double dg;
	if (smax > stol)
		dg = beta;
	else if (smax < -stol)
		dg = gam;
	else {
		ASDCDP3D_COUNT(macaulay_kink);
		dg = 0.5 * (beta + gam);
	}
	// divided and not multiplied by the reciprocal: x/(1-a) and x*(1/(1-a)) are
	// not the same double, and this expression enters n, which enters den, which
	// the whole iteration is a ratio of
	double den = 1.0 - alpha;
	for (int i = 0; i < 3; ++i)
		out(i) = (dq(i) + alpha + dg * P(i)) / den;
	for (int i = 3; i < 6; ++i)
		out(i) = (dq(i) + dg * P(i)) / den;
}

void ASDPlasticDamageConcrete3DMaterial::flowDirection(const Vector& s, Vector& out) const
{
	static Vector dq(6);
	double t = std::tan(psi);
	double q = misesOf(s);
	double a0 = ecc * hard_t.initialYield() * t;
	double t3 = t / 3.0;
	if (!std::isfinite(q)) {
		// misesOf squares the deviator, so it overflows above 1e154 while the
		// RATIO it is wanted for is perfectly well behaved: an unrepresentably
		// large q is by any measure much larger than a0, so q/sqrt(a0^2+q^2) is
		// 1. Dropping the deviatoric term instead - which is what testing the
		// denominator did - silently turns the flow purely volumetric at large
		// stress
		ASDCDP3D_COUNT(mises_overflow);
		dqds(s, dq);
		for (int i = 0; i < 3; ++i)
			out(i) = t3 + dq(i);
		for (int i = 3; i < 6; ++i)
			out(i) = dq(i);
		return;
	}
	double den = std::sqrt(a0 * a0 + q * q);
	if (den > 0.0) {
		dqds(s, dq);
		double f = q / den;
		for (int i = 0; i < 3; ++i)
			out(i) = t3 + f * dq(i);
		for (int i = 3; i < 6; ++i)
			out(i) = f * dq(i);
	}
	else {
		// AT psi = 0 the potential is purely deviatoric, a0 = 0 with it, and on
		// the hydrostatic axis m is identically ZERO: there is then no flow
		// direction at all and no plastic strain can relieve a hydrostatic
		// stress. integrate() detects that as a step with no admissible
		// multiplier rather than dividing by it
		out.Zero();
		for (int i = 0; i < 3; ++i)
			out(i) = t3;
	}
}

double ASDPlasticDamageConcrete3DMaterial::splitWeight(void) const
{
	// ON THE EFFECTIVE STRESS AND NOT ON THE NOMINAL ONE - see the header for
	// what reading it on the nominal one costs
	static Vector d(3);
	if (ASDSpectralSplit::eigenvalues(sbar, d) < 0) {
		ASDCDP3D_COUNT(eigen_error);
		return 0.5;
	}
	return positiveFraction(d);
}

void ASDPlasticDamageConcrete3DMaterial::hardeningRates(const Vector& m, double r,
	double& h_t, double& h_c) const
{
	static Vector d(3);
	if (ASDSpectralSplit::eigenvalues(m, d) < 0) {
		ASDCDP3D_COUNT(eigen_error);
		h_t = h_c = 0.0;
		return;
	}
	// the eigenvalues arrive DESCENDING, so d(0) is the maximum and d(2) the
	// minimum. The Macaulay brackets are what keeps both measures
	// non-decreasing, which is a thermodynamic requirement and not a guard
	h_t = r * macaulay(d(0));
	h_c = (1.0 - r) * macaulay(-d(2));
}

// ==================================================================== //
//  6  the reductions and the damaged elastic law                       //
// ==================================================================== //

double ASDPlasticDamageConcrete3DMaterial::omega(double kappa, double df,
	const ASDCDPHardeningCurve& hard) const
{
	if (df <= 0.0)
		return 1.0;
	double q, dq;
	hard.evaluate(kappa, q, dq);
	double z = df * E * kappa + q;
	if (!(z > 0.0) || !std::isfinite(z)) {
		ASDCDP3D_COUNT(omega_guard);
		return 1.0;
	}
	// BOUNDED IN (0, 1] BY CONSTRUCTION, because df*E*kappa >= 0 - and the clamp
	// on the way out is there so that the BOUND, not the arithmetic, is what the
	// rest of the model relies on: W is the operator every stress and every
	// tangent goes through, and a reduction that came out at 1 + 1e-16 or at
	// -0.0 from a rounded z would put that error in all of them
	double w = q / z;
	if (w < 0.0) w = 0.0;
	if (w > 1.0) w = 1.0;
	return w;
}

double ASDPlasticDamageConcrete3DMaterial::domega(double kappa, double df,
	const ASDCDPHardeningCurve& hard) const
{
	if (df <= 0.0)
		return 0.0;
	double q, dq;
	hard.evaluate(kappa, q, dq);
	double z = df * E * kappa + q;
	double z2 = z * z;
	// z*z is the one squaring in this file that can leave double precision from
	// BELOW: z underflows its own square at 1.5e-154, and a 0.0 there turns a
	// first-order term of the denominator into an inf. The test is therefore on
	// z*z being usable, not on z
	if (!(z > 0.0) || !(z2 > 0.0) || !std::isfinite(z2)) {
		ASDCDP3D_COUNT(omega_guard);
		return 0.0;
	}
	return df * E * (dq * kappa - q) / z2;
}

double ASDPlasticDamageConcrete3DMaterial::leeFenvesReduction(double wt,
	double wc, double r) const
{
	// (1-d) = (1 - s_t*d_c)*(1 - s_c*d_t), Lee and Fenves (1998), and the same
	// expression Abaqus' CDP uses. See DamageCombination in the header for the
	// two limits that define w_t and w_c and for the pair at which this meets
	// the Faria split.
	//
	// NO CLAMP, deliberately, unlike omega(): here the bound is ALGEBRAIC. Both
	// s are 1 minus a product of two numbers in [0,1] and both d are 1 minus an
	// omega that is already clamped into (0,1], so every factor is in [0,1] and
	// the product cannot leave (0,1] by rounding the way a ratio q/z can
	double d_t = 1.0 - wt;
	double d_c = 1.0 - wc;
	double s_t = 1.0 - stiffness_recovery_t * r;
	double s_c = 1.0 - stiffness_recovery_c * (1.0 - r);
	return (1.0 - s_t * d_c) * (1.0 - s_c * d_t);
}

void ASDPlasticDamageConcrete3DMaterial::leeFenvesDerivatives(double wt,
	double wc, double r, double dwt, double dwc, double& g_t, double& g_c,
	double& g_r) const
{
	// The chain:
	//
	//   d(1-d)/dkt = (1 - s_t*d_c) * s_c * domega_t
	//   d(1-d)/dkc = (1 - s_c*d_t) * s_t * domega_c
	//   d(1-d)/dr  = w_t*d_c*(1 - s_c*d_t) - w_c*d_t*(1 - s_t*d_c)
	//
	// because d_t = 1 - omega_t makes d d_t/dkt = -domega_t, and the minus meets
	// the one in front of s_c; and because ds_t/dr = -w_t while ds_c/dr = +w_c,
	// which is why the third one carries the two weights and has no fixed sign.
	//
	// Both domega are negative wherever the secant modulus decreases, and both s
	// are non-negative, so g_t and g_c are negative: the scalar reduction falls
	// as either measure advances, which is what makes the damage terms of the
	// denominator help rather than fight - exactly as in the Faria branch, and
	// for the same reason. g_r is the term the header says was measured rather
	// than assumed away
	double d_t = 1.0 - wt;
	double d_c = 1.0 - wc;
	double s_t = 1.0 - stiffness_recovery_t * r;
	double s_c = 1.0 - stiffness_recovery_c * (1.0 - r);
	g_t = (1.0 - s_t * d_c) * s_c * dwt;
	g_c = (1.0 - s_c * d_t) * s_t * dwc;
	g_r = stiffness_recovery_t * d_c * (1.0 - s_c * d_t)
		- stiffness_recovery_c * d_t * (1.0 - s_t * d_c);
}

double ASDPlasticDamageConcrete3DMaterial::splitWeightRate(const Vector& Cm,
	double pf, double r) const
{
	// r = P/T with P = sum <sbar_i> and T = sum |sbar_i|, so its rate needs the
	// rate of each PRINCIPAL VALUE of the effective stress. Along the frozen
	// update sbar moves by -pf*(C:m) per unit multiplier, and the derivative of
	// a simple eigenvalue in a direction is the direction's own quadratic form -
	// v_i . (dsbar) . v_i - which is EXACT and needs no perturbation of the
	// decomposition. d_split and V_split already hold the one this state paid
	// for.
	//
	// THE FLOOR IS THE ONE positiveFraction USES, and it has to be the same
	// number: below it that routine stops reading the state and returns the
	// constant 0.5, so there r does not depend on lambda at all and its rate is
	// zero rather than a large ratio with a tiny denominator
	double T = std::abs(d_split(0)) + std::abs(d_split(1)) + std::abs(d_split(2));
	if (!(T > 1.0e-14))
		return 0.0;
	double dP = 0.0;
	double dT = 0.0;
	for (int j = 0; j < 3; ++j) {
		double v0 = V_split(0, j);
		double v1 = V_split(1, j);
		double v2 = V_split(2, j);
		// the quadratic form on TENSOR shear components: the off-diagonal pairs
		// are counted twice, which is where the factor 2 comes from
		double c = Cm(0) * v0 * v0 + Cm(1) * v1 * v1 + Cm(2) * v2 * v2
			+ 2.0 * (Cm(3) * v0 * v1 + Cm(4) * v1 * v2 + Cm(5) * v0 * v2);
		double ds = -pf * c;
		if (d_split(j) > 0.0)
			dP += ds;
		dT += (d_split(j) >= 0.0) ? ds : -ds;
	}
	return (dP - r * dT) / T;
}

int ASDPlasticDamageConcrete3DMaterial::effective(const Vector& eps,
	const Vector& ep, double kt_, double kc_, Vector& out)
{
	static Vector de(6);
	static Matrix PT(6, 6);
	static Matrix PC(6, 6);
	double mu2 = E / (1.0 + nu);
	double lam = nu * mu2 / (1.0 - 2.0 * nu);
	// sbar = C : (eps - eps_p), with eps_p the PLASTIC strain. ep_cr, the
	// cracking accumulator, is NOT part of the elastic law at all
	for (int i = 0; i < 6; ++i)
		de(i) = eps(i) - ep(i);
	elasticStress(lam, mu2, de, sbar);
	// the split, and the record it comes from
	if (ASDSpectralSplit::spectral(sbar, d_split, V_split) < 0) {
		ASDCDP3D_COUNT(eigen_error);
		return EC_Eigen_Error;
	}
	// THE BIFURCATION, and the Faria block below is left VERBATIM rather than
	// factored against the other one even where the two look alike. The default
	// branch has to stay bit-identical to what the whole robustness campaign was
	// measured on, and a shared expression that a compiler may re-associate is
	// exactly how that guarantee is lost for no gain
	if (damage_combination == DC_Faria) {
		ASDSpectralSplit::splitFromSpectral(d_split, V_split, PT, PC);
		double wt = omega(kt_, damage_t, hard_t);
		double wc = omega(kc_, damage_c, hard_c);
		// W = omega_t*PT + omega_c*PC
		W.addMatrix(0.0, PT, wt);
		W.addMatrix(1.0, PC, wc);
		sbar_pos.addMatrixVector(0.0, PT, sbar, 1.0);
		sbar_neg.addMatrixVector(0.0, PC, sbar, 1.0);
		// the NOMINAL stress: what equilibrium sees and what the surface is written in
		out.addMatrixVector(0.0, W, sbar, 1.0);
	}
	else {
		// LEE-FENVES: ONE SCALAR, and the nominal stress coaxial with sbar.
		//
		// r IS READ OFF THE DECOMPOSITION JUST PAID FOR, not through
		// splitWeight(): the two are the same number to the bit - both routines
		// of ASDSpectralSplit fill the same matrix from the same Voigt vector
		// and call the same eigen3 - so this is one decomposition instead of
		// two, and it also removes any way for the r that builds the operator to
		// differ from the r the return mapping splits the flow with.
		//
		// NEITHER PROJECTOR IS BUILT AT ALL. PT/PC exist in the Faria branch to
		// carry a DIRECTION, and this branch has none; sbar_pos / sbar_neg are
		// read by corrector() alone, which takes the sbar-wide term here. So
		// this branch is CHEAPER than the default, not dearer: no
		// splitFromSpectral, no two 6x6 assemblies, no two matrix-vector
		// products. They are zeroed rather than left alone so that nothing
		// downstream - the peek record, the serialization - can carry a stale
		// value that belonged to another state
		double r = positiveFraction(d_split);
		double wt = omega(kt_, damage_t, hard_t);
		double wc = omega(kc_, damage_c, hard_c);
		double w = leeFenvesReduction(wt, wc, r);
		W.Zero();
		for (int i = 0; i < 6; ++i)
			W(i, i) = w;
		sbar_pos.Zero();
		sbar_neg.Zero();
		for (int i = 0; i < 6; ++i)
			out(i) = w * sbar(i);
	}
	return 0;
}

void ASDPlasticDamageConcrete3DMaterial::rebuildOperator(const Vector& d_rec,
	const Matrix& V_rec, double kt_, double kc_)
{
	static Matrix PT(6, 6);
	static Matrix PC(6, 6);
	if (damage_combination == DC_Faria) {
		ASDSpectralSplit::splitFromSpectral(d_rec, V_rec, PT, PC);
		double wt = omega(kt_, damage_t, hard_t);
		double wc = omega(kc_, damage_c, hard_c);
		W.addMatrix(0.0, PT, wt);
		W.addMatrix(1.0, PC, wc);
	}
	else {
		// r IS FROZEN WITH THE SPLIT IT COMES FROM, and that is the whole reason
		// this branch reads it off the record instead of off the current state.
		//
		// Under IMPL-EX the projectors are frozen because a spectral
		// decomposition is not a monotone function of an irreversible measure;
		// r is a function of the SAME decomposition, so re-reading it on the
		// extrapolated stress would freeze half of one object and extrapolate
		// the other half. It would also cost the property that pays for the
		// whole scheme: with r frozen, kt and kc are constants of the step, so
		// the scalar reduction is a constant of the step and sigma stays AFFINE
		// in eps - which is what makes W:C the exact algorithmic tangent rather
		// than an approximation of one. Read on the moving stress, r would
		// depend on eps and that statement would be false.
		//
		// It is read here and not taken from r_commit because r_commit is a
		// diagnostic of the committed STEP - zeroed when the step was elastic -
		// while this is a property of the committed STATE
		double r = positiveFraction(d_rec);
		double w = leeFenvesReduction(omega(kt_, damage_t, hard_t),
			omega(kc_, damage_c, hard_c), r);
		W.Zero();
		for (int i = 0; i < 6; ++i)
			W(i, i) = w;
	}
}

// ==================================================================== //
//  7  the return mapping                                               //
// ==================================================================== //

double ASDPlasticDamageConcrete3DMaterial::residualScale(void) const
{
	// THE STRENGTH SCALE OF THE PROBLEM, AND IT MUST NOT MOVE WITH THE STATE.
	//
	// It used to be built from the CURRENT hardening variables,
	//
	//     scale = max(qc, 1e-12) * (1 + |beta(qt, qc)|/(1-alpha))
	//
	// which made tol a relative tolerance against a moving target. Two effects
	// fought in there: qc falling under softening TIGHTENS it, while qt -> 0
	// blows beta up and LOOSENS it, roughly as qc^2/qt. Measured over the
	// 55-case suite, the loosening wins by four orders of magnitude - the
	// absolute tolerance actually used travelled from 7e-9 to 3e-5 while this
	// fixed reference is 3e-9 throughout - and that loose end is the same order
	// as an independent admissibility check's own threshold, which is exactly
	// where the material was accepting states that are outside the surface.
	//
	// Measured, replacing it with the peak strengths: SILENT VIOLATIONS 2 -> 0
	// (a state accepted while inadmissible is the one failure mode this model
	// must not have), failed steps 409 -> 408, so it does not even cost the
	// convergence it was buying. On every case that converges cleanly the answer
	// moves by at most 2.4e-8 MPa, 8e-10 of fc - it changes where the iteration
	// stops, not what it converges to.
	//
	// stressReference() is the peak of both backbones, fixed at construction, and
	// it is what the IMPL-EX error already normalizes against - so the two
	// tolerances of this material are now measured against the same yardstick.
	//
	// NOTE FOR THE PORT: the tolerance decides WHERE the iteration stops, so the
	// Python bench needs the same formula or the bit-agreement that validates
	// this material is gone.
	//
	// tol stays dimensionless either way: verified by re-running the whole suite
	// with every stress scaled by exact powers of two, bit-identical results and
	// identical decisions from fc = 0.029 to fc = 3.1e7.
	return stressReference();
}

double ASDPlasticDamageConcrete3DMaterial::residualStress(double f, double qt,
	double qc) const
{
	// See the header for what this is for and, more importantly, for what it is
	// NOT for: it decides the VERDICT and never the acceptance
	double amp = 1.0 + std::abs(surfaceBeta(qt, qc)) / (1.0 - alpha);
	return (amp > 0.0) ? std::abs(f) / amp : std::abs(f);
}

void ASDPlasticDamageConcrete3DMaterial::corrector(double r, const Vector& m,
	double h_t, double h_c, double kt_, double kc_, Vector& out) const
{
	static Vector Cm(6);
	double mu2 = E / (1.0 + nu);
	double lam = nu * mu2 / (1.0 - 2.0 * nu);
	elasticStress(lam, mu2, m, Cm);
	// the plastic term: W : C : m, scaled by how much of the flow is permanent.
	// The scaling is applied AFTER the product and not folded into
	// addMatrixVector's factor, which would multiply each component of Cm before
	// the sum and give a different last bit
	out.addMatrixVector(0.0, W, Cm, 1.0);
	double pf = plasticShare(r);
	for (int i = 0; i < 6; ++i)
		out(i) *= pf;
	// AND THE TWO DAMAGE TERMS. Both domega are negative, so both ADD. They are
	// what lets the return mapping bring a state back without moving any strain,
	// which is the whole content of df > 0; at df = 0 the corresponding domega is
	// identically zero and the term disappears.
	//
	// THIS IS THE PART THE COMBINATION CHANGES, and it is the only part of the
	// corrector that does: the plastic term above is pf*(W:C:m) in both branches,
	// because W is the operator either way and the flow is split by r either way.
	// What differs is WHERE the damage terms push. In the Faria branch each
	// reduction moves its own spectral part, so the two terms ride on sbar_pos
	// and sbar_neg; in the Lee-Fenves branch there is one scalar in front of the
	// whole stress, so both terms ride on sbar itself and the return they produce
	// is RADIAL by construction
	if (damage_combination == DC_Faria) {
		double a_t = domega(kt_, damage_t, hard_t) * h_t;
		double a_c = domega(kc_, damage_c, hard_c) * h_c;
		for (int i = 0; i < 6; ++i)
			out(i) -= a_t * sbar_pos(i);
		for (int i = 0; i < 6; ++i)
			out(i) -= a_c * sbar_neg(i);
	}
	else {
		double g_t, g_c, g_r;
		leeFenvesDerivatives(omega(kt_, damage_t, hard_t),
			omega(kc_, damage_c, hard_c), r,
			domega(kt_, damage_t, hard_t), domega(kc_, damage_c, hard_c),
			g_t, g_c, g_r);
		// AND THE THIRD TERM, r's own rate, which the Faria branch has no
		// counterpart for because r is not in its operator. Cm is already here
		// and is exactly the direction sbar moves along, so the whole term costs
		// three quadratic forms
		double a = g_t * h_t + g_c * h_c + g_r * splitWeightRate(Cm, pf, r);
		for (int i = 0; i < 6; ++i)
			out(i) -= a * sbar(i);
	}
}

double ASDPlasticDamageConcrete3DMaterial::lambdaCap(double h_t, double h_c,
	const Vector& m, const Vector& eps) const
{
	// neither hardening measure may leave the ten-fold neighbourhood of the
	// curve that indexes it: the measures are bounded by the curves they index,
	// and past that the iteration is diverging and not converging slowly
	double kt_bound = 10.0 * hard_t.lastKappa() + 1.0;
	double kc_bound = 10.0 * hard_c.lastKappa() + 1.0;
	double cap = std::numeric_limits<double>::infinity();
	if (h_t > 0.0) {
		cap = std::min(cap, (kt_bound - kt_commit) / h_t);
		ASDCDP3D_COUNT(cap_binds_kt);
	}
	if (h_c > 0.0) {
		cap = std::min(cap, (kc_bound - kc_commit) / h_c);
		ASDCDP3D_COUNT(cap_binds_kc);
	}
	if (!std::isfinite(cap)) {
		// a step where neither measure advances, so the statement above is
		// vacuous: the plastic strain may not exceed ten times the total strain,
		// which is the only other scale in the problem
		ASDCDP3D_COUNT(cap_from_strain);
		double mn = maxAbs(m);
		if (mn <= 0.0) {
			ASDCDP3D_COUNT(cap_zero);
			return 0.0;    // no flow direction at all: nothing to try
		}
		cap = 10.0 * std::max(maxAbs(eps), 1.0e-6) / mn;
	}
	return std::max(cap, 0.0);
}

bool ASDPlasticDamageConcrete3DMaterial::atApex(const Vector& s) const
{
	// THE VERTEX IS WHERE THE DEVIATOR HAS NOTHING LEFT TO SAY: mean stress
	// tensile, and the Mises equivalent no larger than it. Both apex populations
	// satisfy that with room to spare (hydrostatic tension has q = 0
	// identically) and every harmful den <= 0 event measured on the reference
	// implementation misses it by at least a factor of 24, so the test separates
	// them without being delicate
	double p = traceOf(s) / 3.0;
	return p > 0.0 && misesOf(s) <= p;
}

void ASDPlasticDamageConcrete3DMaterial::trialAt(double lam, const Vector& eps,
	const Vector& m, double pf, double h_t, double h_c, Trial& out)
{
	double f_cr = (1.0 - pf) * lam;
	double f_pl = pf * lam;
	for (int i = 0; i < 6; ++i) {
		out.ep_cr(i) = ep_cr_commit(i) + f_cr * m(i);
		out.ep_pl(i) = ep_pl_commit(i) + f_pl * m(i);
	}
	// the two accumulators always sum to lambda*m: what is not permanent is
	// carried by the reductions
	out.kt = std::max(kt_commit, kt_commit + lam * h_t);
	out.kc = std::max(kc_commit, kc_commit + lam * h_c);
	out.ok = (effective(eps, out.ep_pl, out.kt, out.kc, out.s) == 0);
	double q1, dq1, q2, dq2;
	hard_t.evaluate(out.kt, q1, dq1);
	hard_c.evaluate(out.kc, q2, dq2);
	++n_residual;
	ASDCDP3D_COUNT(residuals);
	out.f = out.ok ? yieldFunction(out.s, q1, q2)
		: std::numeric_limits<double>::quiet_NaN();
	out.ftol = tol * residualScale();
	if (!std::isfinite(out.f))
		out.ok = false;
}

bool ASDPlasticDamageConcrete3DMaterial::bracket(double f0, double lam0, double cap,
	const Vector& eps, const Vector& m, double pf, double h_t, double h_c,
	double& lam_out, Trial& out)
{
	ASDCDP3D_COUNT(bracket_called);
	// THE BEST ITERATE THE SCAN SEES, kept so that a bracket which finds NO root
	// still hands back something better than nothing. The trials are paid for
	// either way; see the note on the failure path below for what they are worth.
	static Trial keep;
	double keep_lam = lam0;
	double keep_abs = std::numeric_limits<double>::infinity();
	bool keep_any = false;
	if (!(cap > lam0))
		return false;
	// THE SCAN IS GEOMETRIC AND UPWARD, because the multiplier a step needs
	// spans the whole range: the elastic-limit steps of a fine ramp want 1e-9
	// and an apex event wants the largest admissible value. A linear scan would
	// need a million points to see the first and a geometric one sees both in
	// twenty
	double lo = lam0;
	double flo = f0;
	bool found = false;
	double hi = cap;
	for (int k = 20; k >= 0; --k) {
		double lam = lam0 + (cap - lam0) * std::ldexp(1.0, -k);
		trialAt(lam, eps, m, pf, h_t, h_c, out);
		if (!out.ok)
			break;
		if (std::abs(out.f) < keep_abs) {
			keep_abs = std::abs(out.f);
			keep_lam = lam;
			keep = out;
			keep_any = true;
		}
		if (std::abs(out.f) <= out.ftol) {
			ASDCDP3D_COUNT(bracket_scan_hit);
			ASDCDP3D_COUNT(bracket_found);
			lam_out = lam;
			return true;                       // landed on it outright
		}
		if (out.f * flo < 0.0) {
			hi = lam;
			found = true;
			break;
		}
		lo = lam;
		flo = out.f;
	}
	if (!found) {
		// NO SIGN CHANGE: there is no admissible multiplier along this direction
		// and saying so is the right answer. But the scan WALKED PAST states far
		// better than the one the caller will otherwise commit, and throwing
		// them away is not part of that answer.
		//
		// Measured on the worst step of the Lee-Fenves branch, 'psi = 0' step 0
		// from a virgin material: F along the flow direction reads 83.30, 84.76,
		// 67.49, 2.78, 74183 - never negative, so no root - and the step used to
		// commit the elastic predictor at 83.3 while the scan had stood at 2.78.
		// Over that branch's 55 cases the worst radial overshoot goes from 20.58
		// to 5.04 MPa, and the DEFAULT branch does not move one bit: there the
		// bracket almost always finds its root, and when it does the adopted
		// state is inside ftol, which no scan point can beat.
		if (keep_any) {
			lam_out = keep_lam;
			out = keep;
		}
		return false;
	}
	ASDCDP3D_COUNT(bracket_found);
	for (int it = 0; it < 30; ++it) {
		double mid = 0.5 * (lo + hi);
		trialAt(mid, eps, m, pf, h_t, h_c, out);
		if (!out.ok)
			return false;
		if (std::abs(out.f) <= out.ftol) {
			lam_out = mid;
			return true;
		}
		if (out.f * f0 > 0.0)
			lo = mid;
		else
			hi = mid;
	}
	// the interval is exhausted: hand back the end that is INSIDE the surface,
	// so a step that ran out of bisections still returns an admissible state
	// rather than one that overshot
	ASDCDP3D_COUNT(bracket_exhausted);
	lam_out = (f0 > 0.0) ? hi : lo;
	trialAt(lam_out, eps, m, pf, h_t, h_c, out);
	return out.ok;
}

int ASDPlasticDamageConcrete3DMaterial::integrate(const Vector& eps)
{
	// Cutting plane, with the damage terms in the denominator.
	//
	// THE UPDATE IS BACKWARD EULER, NOT ACCUMULATED, and that is not a
	// refinement - it is what makes the model work at all. A textbook cutting
	// plane adds dlambda_i * m_i at every iteration, which makes every measure a
	// QUADRATURE OF THE ITERATION PATH. In uniaxial compression with dilatancy
	// that path is not uniaxial: the elastic predictor sees the lateral
	// expansion the plastic flow is about to produce, WITHOUT the flow, so its
	// transverse effective stress is POSITIVE - measured 0.39 to 1.23 MPa
	// against an axial 17.8 to 24.3, on 99 of 271 iterates. Accumulating
	// through those states drove the cracking measure to kt = 1.5e-4 on a path
	// that never left compression, took omega_t to 0.32, and put the zero-stress
	// point of the compressive unloading at -4.7e-4 instead of the plastic
	// strain -2.03e-3. Writing the update with the TOTAL multiplier and the
	// direction at the CURRENT iterate removes the quadrature: at convergence
	// the rates are those of the converged state, where r = 0 and the tensile
	// measure does not move.
	//
	// THE CORRECTOR IS DAMPED AND THE FALLBACK IS A BRACKET, and those two are
	// what make the difference between an algorithm and a formula. den is NOT
	// wrong - it is the derivative of the update with m, r, pf, h_t and h_c
	// FROZEN at the current iterate, which is what a cutting plane is, and
	// against a central difference of precisely that frozen function it agrees
	// to 0.2%. But nothing makes the frozen quantities mild functions of the
	// state, and along a stalling step r went 0.998 -> 0.573 -> 0.093 -> 0.0007.
	// Once the iterate leaves the basin where they are nearly constant the
	// frozen derivative stops describing the function being solved: measured,
	// den stayed at +1e4 while the true -dF/dlambda turned -554, and the
	// corrector then walked UPHILL for thirty consecutive iterations, |F| growing
	// from 5e-5 to 2e4. den <= 0 cannot catch that - it tests the frozen
	// Jacobian, which is healthy and positive throughout. What catches it is
	// testing the thing one actually wants: THE RESIDUAL HAS TO GO DOWN.
	//
	// (That particular runaway was a BUG and is gone: r was being read off the
	// damaged stress. With it corrected the globalization is never exercised on
	// the tension-tail, alternating and confined-shear protocols - all three
	// return bit-identical results with the line search disabled, because the
	// full Newton step is accepted every time. It is still load-bearing for
	// exactly one state, the apex, where the shipped loop takes 5 failures
	// against the undamped 16.)
	ASDCDP3D_COUNT(integrate);

	static Vector n(6);
	static Vector m(6);
	static Vector A(6);
	static Trial st;
	static Trial adopted;
	// the best iterate seen, restored if the loop gives up: a failed step must
	// not poison the state it failed from. Without this an apex event left
	// kc = 6.59 on a curve that ends at 1.4e-2 and every later step inherited
	// it. The elastic predictor is IN the ranking, at lambda = 0, and it wins
	// only when no trial was ever accepted - which is the honest outcome of a
	// step where no multiplier reduces the residual
	static Vector best_ep_cr(6);
	static Vector best_ep_pl(6);
	double best_f, best_kt, best_kc, best_lam;

	ep_cr = ep_cr_commit;
	ep_pl = ep_pl_commit;
	kt = kt_commit;
	kc = kc_commit;
	n_iter = 0;
	dlambda = 0.0;
	plastic = false;
	failed = false;
	stagnated = false;
	bisected = 0;
	n_residual = 1;
	if (effective(eps, ep_pl, kt, kc, stress) != 0)
		return EC_Eigen_Error;
	double qt, dqt, qc, dqc;
	hard_t.evaluate(kt, qt, dqt);
	hard_c.evaluate(kc, qc, dqc);
	double f = yieldFunction(stress, qt, qc);
	if (!std::isfinite(f))
		return EC_Eigen_Error;
	if (f <= tol * residualScale()) {
		ASDCDP3D_COUNT(elastic);
		return 0;
	}

	ASDCDP3D_COUNT(plastic);
	plastic = true;
	best_f = std::abs(f);
	best_ep_cr = ep_cr;
	best_ep_pl = ep_pl;
	best_kt = kt;
	best_kc = kc;
	best_lam = 0.0;
	// THE WATCHDOG'S REFERENCE, and it has to be the residual this step STARTED
	// from rather than the best one seen since. Referencing the running best
	// makes the bound tighten as the iteration succeeds, which sounds right and
	// forbids the one thing the cutting plane is documented to need: a corrector
	// that overshoots and is walked back by a NEGATIVE one, which happens
	// whenever a hardening modulus changes inside the step. At a fine
	// discretization the entry residual is small, so the running best is small,
	// so a single legitimate overshoot breaks a bound relative to it - measured
	// on uniaxial compression at n = 400, the backbone error went to 2.0e-4 with
	// 12 failed steps against 5.6e-8 and none
	double f_entry = std::abs(f);
	bool broke = false;
	// dlambda is parameterizing kappa_t and not the flow: see the no-flow branch
	bool on_damage = false;

	for (int it = 1; it <= max_iter; ++it) {
		n_iter = it;
		ASDCDP3D_COUNT(iterations);
		double smax = 0.0;
		dFds(stress, qt, qc, n, smax);
		flowDirection(stress, m);
		double r = splitWeight();
		double h_t, h_c;
		hardeningRates(m, r, h_t, h_c);
		hard_t.evaluate(kt, qt, dqt);
		hard_c.evaluate(kc, qc, dqc);
		double dF_dkt, dF_dkc;
		dFdKappa(smax, qt, qc, dqt, dqc, dF_dkt, dF_dkc);
		double pf = plasticShare(r);
		corrector(r, m, h_t, h_c, kt, kc, A);
		double den = contract(n, A) - (dF_dkt * h_t + dF_dkc * h_c);
		double cap = lambdaCap(h_t, h_c, m, eps);
		double lam0 = dlambda;
		bool got = false;
		double got_lam = 0.0;
		// invalidated so that the ranking below cannot mistake a trial left by a
		// PREVIOUS iteration for one this iteration produced
		adopted.ok = false;

		// ---- THE ITERATION CAME BACK TO WHERE IT STARTED: it is CYCLING ---- //
		//
		// den > 0 the whole way, so nothing below ever runs, and the corrector
		// simply goes round. Measured on 'implex, dials 0.5/0.5' step 9, and the
		// trace is unambiguous - period three, repeated until max_iter:
		//
		//   it=1  f= 16.7229  lam0=0         r=1      h_t=0.214  h_c=0
		//   it=2  f=  4.98326 lam0=0.001717  r=0.937  h_t=0.649  h_c=0.0506
		//   it=3  f=-13.6788  lam0=0.001934  r=1      h_t=0.199  h_c=0
		//   it=4  f= 16.7229  lam0=0         ...  identical to it=1
		//
		// At the third iterate F is NEGATIVE - the state is inside the surface -
		// so Newton walks back, the accumulated multiplier would go negative, the
		// clamp puts it at zero, and that is the starting point again. It cycles
		// because the frozen quantities are not mild functions of the state: the
		// split weight r swings 1 -> 0.937 -> 1 between iterates and takes h_t
		// and h_c with it, so each iterate solves a slightly different problem.
		// best_f ends up at 4.98, which is exactly the residual that step commits.
		//
		// A SECOND ARRIVAL AT lambda == 0 IS THE SIGNATURE and it cannot mean
		// anything else: the first iterate is only accepted with a strictly
		// positive multiplier (f > 0 and den > 0 give dlam > 0), so coming back to
		// zero later means the corrector walked out and was clamped back.
		//
		// Bracketing FROM ZERO is then both available and correct - the scan
		// covers the whole range from below, so the first sign change is the
		// smallest admissible multiplier - and the root is there: on that step F
		// goes +8.49 at lambda 9.3e-4 to -4.01 at 5.3e-3, crossing in
		// [2.63e-3, 3.72e-3].
		//
		// Measured at damageT 1 / damageC 0.3: failed steps 118 -> 92 and the
		// worst F over the suite 1.30 -> 0.698 MPa, 2.3% of fc. At the suite's own
		// dials with damage_t = 0 excluded: 137 -> 116 and 4.98 -> 3.69. The
		// worst deviation on every step the Python bench solves is UNCHANGED on
		// both - same value, same step - so nothing the oracle can vouch for
		// moves. (The median of what still fails rises, 4.3e-5 -> 3.3e-3: the
		// steps this converts were the small-residual ones, and what is left is
		// the harder remainder.)
		bool cycled = (it > 1 && lam0 == 0.0);
		if (cycled && bracket(f, 0.0, cap, eps, m, pf, h_t, h_c, got_lam, adopted)) {
			ASDCDP3D_COUNT(cycle_bracketed);
			got = true;
			++bisected;
		}
		else if (std::isfinite(den) && den > 0.0) {
			ASDCDP3D_COUNT(newton_ok);
			// ---- the Newton corrector, halved until it is acceptable ---- //
			double dlam = f / den;
			// the ACCUMULATED multiplier is what may not go negative - a single
			// corrector may, and does whenever a hardening modulus changes
			// inside the step
			if (lam0 + dlam < 0.0) {
				ASDCDP3D_COUNT(negative_corrector);
				// ---- THE OVERSHOOT HAS A KNOWN INTERVAL: (0, lam0) ---- //
				//
				// f < 0 at lam0 > 0 means the corrector went PAST the root, and
				// the root is then bracketed by two values this step already
				// owns: the residual at lambda = 0 is the elastic predictor's,
				// positive by definition of a plastic step, and the residual
				// here is negative. bracket() is exactly the routine for an
				// interval with a sign change and no usable derivative.
				//
				// WHAT THE CLAMP DOES INSTEAD is throw lambda back to zero, and
				// the limit-cycle rule then brackets on [0, cap] on the NEXT
				// iterate. Two things are wrong with arriving that way: the
				// round trip, and cap is not lam0 - bracket() scans 21 geometric
				// points, so on [0, cap] with cap >> lam0 the resolution near
				// the root Newton just overshot is coarse and the scan can step
				// over a pair of sign changes.
				//
				// WHY IT IS FENCED TO THE LEE-FENVES BRANCH, and it is not
				// conservatism about the algebra - the rule is
				// formulation-agnostic and reads the same either way. It is that
				// the two branches meet this situation at completely different
				// rates: measured over the 55-case suite at damageT 1 /
				// damageC 0.3, negative_corrector fires 2.8 times per 100 plastic
				// steps with the Faria split and 719 with this one, because a
				// scalar reduction collapses the whole stress and F falls far
				// more steeply in lambda. Un-fenced it is not free: the default
				// branch goes from 92 failed steps to 91 with its worst radial
				// overshoot unchanged at 0.1874 MPa, but 255 of its rows stop
				// being bit-identical to the material the whole robustness
				// campaign was measured on, and that reference is worth more
				// than one failed step. Turning it on there is a separate
				// decision with a separate verification, exactly as the
				// normalization flag of ASDSpectralSplit is.
				//
				// MEASURED, at damageT 1 / damageC 0.3, Lee-Fenves: failed steps
				// 900 -> 413 at w = 0/1 and 718 -> 413 at w = 1/1, with the worst
				// radial overshoot of EVERY case unchanged to four figures. So it
				// converts failures without changing the quality of what still
				// fails - which is the direction the other two candidates could
				// not manage: freezing r inside the trial took the worst state
				// from 4.3e3 to 3.5e4 in F, and a second watchdog ceiling on the
				// current residual made BOTH branches worse (the default one 92
				// -> 100/130/165/192 failed steps at G = 10/4/2/1.2).
				if (damage_combination == DC_LeeFenves) {
					static Trial at_zero;
					trialAt(0.0, eps, m, pf, h_t, h_c, at_zero);
					if (at_zero.ok && at_zero.f > 0.0 &&
						bracket(at_zero.f, 0.0, lam0, eps, m, pf, h_t, h_c,
							got_lam, adopted)) {
						ASDCDP3D_COUNT(bracket_below);
						got = true;
						++bisected;
					}
					else {
						// nothing changes sign on the way back: the clamp, as
						// before. The operator is put back on the iterate the
						// line search is about to start from
						trialAt(lam0, eps, m, pf, h_t, h_c, st);
						dlam = -lam0;
					}
				}
				else {
					dlam = -lam0;
				}
			}
			double step = 1.0;
			if (!got)
			for (int bt = 0; bt <= max_backtrack; ++bt) {
				double lam = lam0 + step * dlam;
				if (lam >= 0.0 && lam <= cap) {
					trialAt(lam, eps, m, pf, h_t, h_c, st);
					// A WATCHDOG, NOT A MONOTONE DESCENT TEST, and the gap
					// between those two is where this loop was won or lost. The
					// full Newton step is tried FIRST and accepted unless it is
					// diverging, so a healthy solve takes exactly the corrector
					// the undamped cutting plane took and every tight invariant
					// of the formulation survives bit for bit - kt staying at
					// 1.2e-15 on a purely compressive path, the compressive peak
					// at -fc to 1.2e-10, the tensile backbone to 4.9e-9.
					// Insisting instead on |F| decreasing at EVERY step breaks
					// all three at once (kt to 1.1e-5, the peak to 1.4e-3),
					// because Newton legitimately overshoots the root on its
					// last corrector: that step raises |F| while putting lambda
					// much closer, and the step after it collects the accuracy.
					// Gate it away and the iteration stops one corrector early,
					// every time.
					//
					// What has to be caught is not a single non-monotone step but
					// SUSTAINED divergence, and the third clause is the weakest
					// statement that catches it: a corrector may not leave the
					// state WORSE THAN IT FOUND IT. The runaway crosses that at
					// its fourth step - entry residual 0.215, and |F| passes it
					// on the way from 5e-5 to 2e4 - while an
					// overshoot-and-walk-back stays comfortably inside it
					if (st.ok) {
						double af = std::abs(st.f);
						bool ok_tol = af <= st.ftol;
						bool ok_desc = af < std::abs(f);
						bool ok_dog = af <= ResidualGrowth * f_entry;
						if (ok_tol || ok_desc || ok_dog) {
#ifdef ASDCDP3D_COUNTERS
							if (ok_tol) ASDCDP3D_COUNT(gate_by_tol);
							else if (ok_desc) ASDCDP3D_COUNT(gate_by_descent);
							else ASDCDP3D_COUNT(gate_by_watchdog);
							if (bt > 0) ASDCDP3D_COUNT(backtrack);
#endif
							got = true;
							got_lam = lam;
							adopted = st;
							break;
						}
					}
				}
				else {
					ASDCDP3D_COUNT(lam_out_of_range);
				}
				step *= 0.5;
			}
#ifdef ASDCDP3D_COUNTERS
			if (!got) ASDCDP3D_COUNT(backtrack_exhausted);
#endif
		}
		else {
			ASDCDP3D_COUNT(den_nonpositive);
			// ---- the gradient is unusable AT THE APEX: bracket it ---- //
			//
			// THE APEX GATE STAYS, AND IT IS NOT CONSERVATISM - it was measured.
			// Opening the bracket to every den <= 0 looks like a large win: the
			// failed steps over the 55-case suite go 417 -> 291, 'psi = 0' clears
			// outright, and 51 of 55 cases stay bit-identical. It is still wrong,
			// because bracket() solves F = 0 and F is NOT MONOTONE in lambda:
			// away from the apex its geometric scan can land on a DIFFERENT root
			// from the one Newton was walking to, and the loop then accepts it -
			// admissible, converged, and not the answer. Measured on
			// 'implex, dials 0.5/0.5': at its first step the shipped code agrees
			// with the Python bench to 1.8e-15 while REPORTING the step failed
			// (best_restored happened to hold the converged state), and the
			// ungated bracket converges 12.9 MPa away from it; over that case the
			// gated code is closer to the bench on 150 steps out of 150, worst
			// error 6.3e-4 against 12.9. Trading 163 failed steps for a wrong
			// answer that does not announce itself is the wrong trade.
			//
			// At the apex there is no such ambiguity, which is why the gate is
			// exactly here: the normal is degenerate, the flow has one direction
			// available at most, and there is no second root to be lured onto.
			// ... OR ON THE FIRST ITERATE, WHERE THE SCAN CANNOT PICK THE WRONG
			// ROOT. This is the other half of the story above, and it is what
			// makes the difference between the two.
			//
			// The reason ungating bracket() everywhere was wrong is that it scans
			// STRICTLY UPWARD from lam0: once Newton has overshot, the root it
			// wanted lies BELOW lam0, is invisible to the scan, and the next one
			// up gets adopted instead. At lam0 == 0 that cannot happen - there is
			// nothing below zero - so the first sign change the scan meets IS the
			// smallest admissible multiplier, which is what a return map owes its
			// caller.
			//
			// The step that forced this is the worst one in the whole suite at the
			// dials that will actually be used (damageT 1, damageC 0.3):
			// 'psi = 0' step 0, the first increment of a random walk, where the
			// corrector gives up after ONE iteration with den = -61507 and
			// kt = kc = 0, and commits the elastic predictor at F = 83.3 MPa -
			// 2.8 times fc outside its own yield surface. Scanning F along that
			// same flow direction shows it CHANGES SIGN TWICE well inside the cap
			// (lambda ~ 1.2e-3 and ~7e-3, cap 2.13): the root was there all along
			// and the apex gate was the only thing keeping the one routine that
			// could find it from running. With this branch the step converges to
			// F = 2.7e-9, inside the tolerance.
			//
			// MEASURED, and against the right yardstick. On the 8078 steps the
			// Python bench actually solves, this changes NOTHING - worst
			// deviation 13.56 MPa on both builds, same step - so every state the
			// oracle can vouch for is untouched, and only steps where the bench
			// itself gave up move. At damageT 1 / damageC 0.3: failed steps
			// 133 -> 131 and the worst F over the whole suite 83.3 -> 1.42 MPa,
			// with nothing left above fc. At the suite's own mixed dials:
			// 408 -> 281 failures, median F on what still fails 50.4 -> 0.45 MPa,
			// steps worse than fc 209 -> 79.
			if (atApex(stress) || lam0 == 0.0) {
				// ---- the gradient is unusable AT THE APEX: bracket it ---- //
				ASDCDP3D_COUNT(at_apex);
				if (bracket(f, lam0, cap, eps, m, pf, h_t, h_c, got_lam, adopted)) {
					got = true;
					++bisected;
				}
				// THE TRIGGER IS 'LAMBDA MOVES NOTHING', NOT 'm IS ZERO', and the
			// difference is 13 failed steps at the dials that will be used.
			//
			// m = 0 is only the most obvious way for the multiplier to be inert.
			// The other one has a healthy flow direction and is just as stuck:
			// with damage_t = 1 on a fully tensile state, plasticShare(r = 1) is
			// exactly 0, so lambda creates NO plastic strain; and when both
			// hardening rates are zero it does not advance the measures either.
			// The strain is then unchanged, kappa is unchanged, and F is
			// CONSTANT in lambda - measured on 'hydrostatic tension, psi = 0'
			// step 23, where F reads 1.41991 at every lambda from 7.7e-13 to
			// 0.027, with |m| = 0.866. No bracket can find a root there because
			// there is no dependence to find one in, and - this is the part worth
			// keeping in mind - the increment at that step is 3e-6, so a smaller
			// step cannot help either. It is not a step-size failure.
			//
			// The reduction can still move it, which is what this branch is for.
			// Measured, widening the trigger from 'm = 0' to 'nothing moves':
			// failed steps 131 -> 118 at damageT 1 / damageC 0.3 and the median F
			// of what still fails 3.3e-3 -> 4.3e-5 MPa, with the worst deviation
			// on every step the Python bench solves UNCHANGED - same 13.56 MPa on
			// the same step - so again nothing the oracle can vouch for is touched.
			else if ((h_t == 0.0 && h_c == 0.0
				&& (pf == 0.0 || maxAbs(m) <= 0.0))
				&& (dlambda == 0.0 || on_damage)) {
					// ---- NO FLOW DIRECTION AT ALL: return on the DAMAGE ---- //
					//
					// At psi = 0 the potential is purely deviatoric, and on the
					// hydrostatic axis the deviator is identically zero, so
					// flowDirection returns exactly zero. Then lambda moves NOTHING:
					// not the strain, and not the measures either, because both
					// hardening rates are read off m. There is no multiplier to
					// bracket and no amount of iterating invents one - measured, that
					// is 44 of the 80 failed steps of hydrostatic tension at psi = 0,
					// where bracket() was called 91 times and found a root 12.
					//
					// What CAN relieve it is the reduction. The surface is written in
					// the NOMINAL stress, sigma = omega_t*PT:sbar + omega_c*PC:sbar,
					// so advancing kappa_t at frozen strain scales sigma down; and at
					// the apex sbar is entirely tensile, so that scaling is RADIAL -
					// the one direction the vertex's cone of normals certainly
					// contains, and the only one available.
					//
					// IT NEEDS NO NEW SOLVER. trialAt with m = 0, h_t = 1 and h_c = 0
					// is exactly kappa_t = kappa_t_commit + lambda at frozen strain,
					// so lambda becomes the kappa increment and bracket() drives it
					// unchanged. Nor does it need anything from IMPL-EX: commitState
					// derives its frozen rates by DIVIDING the increments by dlambda,
					// so re-parameterizing is self-consistent - it records
					// kt_rate = 1 with every strain rate zero, which is precisely
					// what this step did.
					//
					// AND THE BRACKET CANNOT BE FOOLED HERE: F falls monotonically
					// towards -qc as omega_t goes to zero, so a sign change exists
					// whenever the reduction can move at all. When it cannot -
					// damage_t = 0, where domega is identically zero - no root is
					// found and the step still reports its honest failure. That is
					// the right answer and not a gap: with no flow AND no reduction
					// the model has no mechanism for that state, and the caller has
					// to hear it rather than be handed a number.
					//
					// The guard is that no plastic strain has moved yet in this step
					// (or that we are already on this branch, where none has):
					// trialAt rebuilds the accumulators from the COMMITTED ones, so
					// re-parameterizing after real plastic flow would silently
					// discard it.
					ASDCDP3D_COUNT(no_flow);
					ASDCDP3D_COUNT(damage_return);
					static Vector no_flow_dir(6);
					no_flow_dir.Zero();
					double cap_k = lambdaCap(1.0, 0.0, no_flow_dir, eps);
					if (bracket(f, lam0, cap_k, eps, no_flow_dir, pf, 1.0, 0.0,
						got_lam, adopted)) {
						ASDCDP3D_COUNT(damage_return_ok);
						got = true;
						on_damage = true;
						++bisected;
					}
				}
			}
		}

		if (!got) {
			// A BRACKET THAT FOUND NO ROOT STILL PAID FOR ITS TRIALS, and it now
			// hands back the best of them - see the failure path of bracket().
			// Ranking it here cannot move a step that succeeds: a step that
			// succeeds adopts a state inside ftol, and no scan point beats that
			if (adopted.ok && std::abs(adopted.f) < best_f) {
				best_f = std::abs(adopted.f);
				best_ep_cr = adopted.ep_cr;
				best_ep_pl = adopted.ep_pl;
				best_kt = adopted.kt;
				best_kc = adopted.kc;
				best_lam = got_lam;
			}
			// NOTHING along this direction reduces the residual. Either the apex
			// with no flow to relieve it, or a kink where the residual has
			// already reached the smallest value it can take. 'best' is what
			// gets restored, so the distinction the flags draw is only about what
			// the caller is told
			// IN STRESS UNITS, not in F units: |F| is a distance times a
			// gradient of ~qc/qt, so comparing it to a strength scale asks a
			// question with the wrong dimensions. See residualStress()
			stagnated = residualStress(best_f, qt, qc)
				<= stagnation_tol * residualScale();
			failed = !stagnated;
#ifdef ASDCDP3D_COUNTERS
			if (stagnated) ASDCDP3D_COUNT(stagnated);
#endif
			broke = true;
			break;
		}

		dlambda = got_lam;
		stress = adopted.s;
		ep_cr = adopted.ep_cr;
		ep_pl = adopted.ep_pl;
		kt = adopted.kt;
		kc = adopted.kc;
		f = adopted.f;
		double ftol = adopted.ftol;
		hard_t.evaluate(kt, qt, dqt);
		hard_c.evaluate(kc, qc, dqc);
		if (std::abs(f) < best_f) {
			best_f = std::abs(f);
			best_ep_cr = ep_cr;
			best_ep_pl = ep_pl;
			best_kt = kt;
			best_kc = kc;
			best_lam = dlambda;
		}
		if (std::abs(f) <= ftol) {
			broke = true;
			break;
		}
	}
	if (!broke) {
		// EXHAUSTING max_iter USED TO BE A FAILURE UNCONDITIONALLY, on the
		// argument that running out of budget says nothing about the problem and
		// there is no evidence the next iteration would not have helped. That
		// argument was right about the reasoning and wrong about the evidence,
		// and the trace of the one step at (1, 0.3) that reaches the cap is what
		// changed it: by the twentieth pass the frozen rates have settled and F
		// sits on a fixed value, so the next iteration demonstrably would NOT
		// have helped - and the best iterate seen was 2.3e-5 MPa from the
		// surface. Calling that a failure costs the host a step cut it does not
		// need.
		//
		// So the cap now asks the SAME question the give-up path asks, and it is
		// the honest one: not "did the loop finish" but "is the state about to be
		// returned admissible". The measure is a stress - see residualStress() -
		// because F is not one.
		ASDCDP3D_COUNT(max_iter);
		stagnated = residualStress(best_f, qt, qc)
			<= stagnation_tol * residualScale();
		failed = !stagnated;
#ifdef ASDCDP3D_COUNTERS
		if (stagnated) ASDCDP3D_COUNT(stagnated);
#endif
	}
	if (failed || stagnated) {
		ASDCDP3D_COUNT(best_restored);
		ep_cr = best_ep_cr;
		ep_pl = best_ep_pl;
		kt = best_kt;
		kc = best_kc;
		dlambda = best_lam;
		if (effective(eps, ep_pl, kt, kc, stress) != 0)
			return EC_Eigen_Error;
	}
	return 0;
}

// ==================================================================== //
//  8  IMPL-EX                                                          //
// ==================================================================== //

double ASDPlasticDamageConcrete3DMaterial::timeFactor(void) const
{
	if (dtime_n_commit <= 0.0)
		return 1.0;
	double f = dtime_n / dtime_n_commit * implex_alpha;
	if (!(f > 0.0)) {
		// CLAMPED AT ZERO. A negative ratio extrapolates an irreversible process
		// BACKWARDS - the plastic strain stops being monotone along the flow
		// direction, the hardening measures decrease and the damage undoes
		// itself - which is not a state this material can be in. Clamping
		// degrades the step to a purely elastic prediction: wrong by O(dt) like
		// everything else here, but a state that EXISTS
#ifdef ASDCDP3D_COUNTERS
		if (f != 0.0) ASDCDP3D_COUNT(time_factor_clamped);
#endif
		return 0.0;
	}
	return f;
}

void ASDPlasticDamageConcrete3DMaterial::extrapolate(const Vector& eps)
{
	ASDCDP3D_COUNT(extrapolate);
	static Vector de(6);
	double dlam = timeFactor() * dlambda_commit;
	for (int i = 0; i < 6; ++i) {
		ep_cr(i) = ep_cr_commit(i) + dlam * mcr_commit(i);
		ep_pl(i) = ep_pl_commit(i) + dlam * m_commit(i);
	}
	kt = kt_commit + dlam * kt_rate_commit;
	kc = kc_commit + dlam * kc_rate_commit;
	// THE SPLIT IS FROZEN AND THE DAMAGE IS EXTRAPOLATED, and the distinction is
	// the whole content of this method. PT/PC are not monotone functions of an
	// irreversible measure - they are a spectral decomposition, they can rotate
	// and exchange eigenvalues freely, and an extrapolated projector means
	// nothing. The two reductions are the opposite: omega is a function of
	// kappa, which IS the monotone measure this scheme extrapolates.
	//
	// AND IT COSTS NOTHING STRUCTURALLY. Freezing the reductions too - which an
	// earlier revision did, on the grounds that sigma must stay AFFINE in eps -
	// buys nothing, because dlam is extrapolated from the COMMITTED state, so kt
	// and kc are CONSTANTS OF THE STEP and neither does omega of them. sigma is
	// exactly as affine either way, and the algorithmic tangent is W:C to the
	// bit. What it buys is the d omega term of the error decomposition, 17 to
	// 23% of the per-step error: the delivered stress error drops about 9% at
	// every refinement, and the COMMITTED trajectory does not move at all
	// because the commit re-solves implicitly anyway
	rebuildOperator(d_split_commit, V_split_commit, kt, kc);
	// the effective stress is kept alongside the nominal one because
	// splitWeight() answers from it and this path never calls effective():
	// without it the IMPL-EX record in commitState() would read an r belonging
	// to whatever state came before
	double mu2 = E / (1.0 + nu);
	double lam = nu * mu2 / (1.0 - 2.0 * nu);
	for (int i = 0; i < 6; ++i)
		de(i) = eps(i) - ep_pl(i);
	elasticStress(lam, mu2, de, sbar);
	stress.addMatrixVector(0.0, W, sbar, 1.0);
	// sbar_pos / sbar_neg are deliberately NOT built here: nothing on the
	// explicit path reads them (corrector() lives inside the return mapping and
	// getTangent() takes the W:C branch under IMPL-EX), and building them would
	// cost two matrix-vector products per step for an output nobody asks for
	dlambda = dlam;
	plastic = dlam > 0.0;
	failed = false;
	stagnated = false;
	bisected = 0;
	n_iter = 0;
	n_residual = 0;
	extrapolated = true;
}

double ASDPlasticDamageConcrete3DMaterial::stressReference(void) const
{
	return std::max(std::max(hard_c.maxEffectiveStress(),
		hard_t.maxEffectiveStress()), 1.0);
}

double ASDPlasticDamageConcrete3DMaterial::implexStressGap(
	const Vector& delivered, const Vector& stress_implicit) const
{
	// THE NORM is the largest Voigt component - the 1D twin takes an absolute
	// value - but THE DENOMINATOR is the same in both, the largest stress the
	// material can carry, or one tolerance would mean different things in 1D and
	// in 3D and a stepper would be comparing numbers that are not comparable
	double gap = 0.0;
	for (int i = 0; i < 6; ++i)
		gap = std::max(gap, std::abs(delivered(i) - stress_implicit(i)));
	return gap / stressReference();
}

double ASDPlasticDamageConcrete3DMaterial::computeImplexErrorMetric(void)
{
	// no extrapolation, no error. This is what makes it safe for the aggregation
	// to call it on everything it has
	if (!implex)
		return 0.0;
	ASDCDP3D_COUNT(peek);
	// the current state holds the EXPLICIT answer: the stress this step
	// delivered to the element, and the state the recorders must keep seeing.
	// stress_implex already holds it - setTrialStrain records it right after the
	// extrapolation - and integrate() below does not write that member, so it
	// survives the throw-away solve without being part of the peek record
	static TrialState delivered;
	saveTrialState(delivered);
	// the implicit answer, at the same trial strain
	if (integrate(strain) < 0) {
		restoreTrialState(delivered);
		// no metric is not a small metric: whoever reads this must not accept
		// the step
		return std::numeric_limits<double>::quiet_NaN();
	}
	double err = implexStressGap(stress_implex, stress);
	// undo. Measuring is not allowed to move the state: the step may still be
	// rejected, and revertToLastCommit() would not put the frozen split back
	restoreTrialState(delivered);
	implex_error = err;
	return err;
}

double ASDPlasticDamageConcrete3DMaterial::implexTimeRatio(void) const
{
	return dtime_0 > 0.0 ? dtime_n / dtime_0 : 1.0;
}

void ASDPlasticDamageConcrete3DMaterial::saveTrialState(TrialState& x) const
{
	x.stress = stress;
	x.sbar = sbar;
	x.sbar_pos = sbar_pos;
	x.sbar_neg = sbar_neg;
	x.ep_cr = ep_cr;
	x.ep_pl = ep_pl;
	x.d_split = d_split;
	x.V_split = V_split;
	x.W = W;
	x.C = C;
	x.kt = kt;
	x.kc = kc;
	x.dlambda = dlambda;
	x.implex_error = implex_error;
	x.n_iter = n_iter;
	x.n_residual = n_residual;
	x.bisected = bisected;
	x.plastic = plastic;
	x.failed = failed;
	x.stagnated = stagnated;
	x.extrapolated = extrapolated;
}

void ASDPlasticDamageConcrete3DMaterial::restoreTrialState(const TrialState& x)
{
	stress = x.stress;
	sbar = x.sbar;
	sbar_pos = x.sbar_pos;
	sbar_neg = x.sbar_neg;
	ep_cr = x.ep_cr;
	ep_pl = x.ep_pl;
	d_split = x.d_split;
	V_split = x.V_split;
	W = x.W;
	C = x.C;
	kt = x.kt;
	kc = x.kc;
	dlambda = x.dlambda;
	implex_error = x.implex_error;
	n_iter = x.n_iter;
	n_residual = x.n_residual;
	bisected = x.bisected;
	plastic = x.plastic;
	failed = x.failed;
	stagnated = x.stagnated;
	extrapolated = x.extrapolated;
}

// ==================================================================== //
//  9  state handling                                                   //
// ==================================================================== //

int ASDPlasticDamageConcrete3DMaterial::setTrialStrain(const Vector& v)
{
	int retval = 0;

	// this material point takes part in the current step, so it takes part in
	// the IMPL-EX error aggregate (see IMPLEXManager.h)
	implexTouch();

	if (!prepare())
		return EC_Generic;

	// save dT
	if (!dtime_is_user_defined) {
		dtime_n = ops_Dt;
		if (!commit_done) {
			dtime_0 = dtime_n;
			dtime_n_commit = dtime_n;
		}
	}

	// ENGINEERING IN, TENSOR INSIDE: one of the only two places the distinction
	// exists - see the header
	strain(0) = v(0);
	strain(1) = v(1);
	strain(2) = v(2);
	strain(3) = 0.5 * v(3);
	strain(4) = 0.5 * v(4);
	strain(5) = 0.5 * v(5);

	if (implex) {
		extrapolate(strain);
		// RECORD WHAT THIS STEP DELIVERS, here and not inside extrapolate(),
		// because commitState calls integrate() over the same members and would
		// erase it: after the commit the only trace of the extrapolated answer
		// would be the scalar error
		stress_implex = stress;
		if (implex_control) {
			// THE LEGACY in-material measurement, read through the
			// NON-DESTRUCTIVE peek and not measured inline, so that the two
			// places that want the error cannot drift. The implicit solve it
			// costs is thrown away: what the step CARRIES is the explicit
			// answer, and what gets committed is a fresh implicit solve at
			// commit time
			double err = computeImplexErrorMetric();
			// and only if the user asked for the old behaviour, fail here
			if (implex_abort_on_error && !(err <= implex_error_tolerance)) {
				if (dtime_n >= implex_time_redution_limit * dtime_0)
					retval = EC_IMPLEX_Error_Control;
			}
		}
	}
	else {
		extrapolated = false;
		retval = integrate(strain);
		// A RETURN MAPPING THAT DID NOT CONVERGE IS REPORTED, and that is not
		// the IMPL-EX policy question. The rule that this material never aborts
		// is about the EXTRAPOLATION ERROR, which is a discretization statement
		// only the analysis can price; a return mapping that ran out of
		// iterations is a local failure with no admissible answer to hand back,
		// and saying so is what lets the integrator cut the step. 'stagnated' is
		// NOT such a case: there the residual stopped at a value that is small
		// on the scale that matters, and the state is accepted
		if (retval == 0 && failed)
			retval = EC_Generic;
		// no extrapolation, so what was delivered IS the implicit answer, and the
		// response means the same thing in both regimes
		stress_implex = stress;
	}

	computeTangent();
	return retval;
}

int ASDPlasticDamageConcrete3DMaterial::setTrialStrain(const Vector& v, const Vector& /*r*/)
{
	return setTrialStrain(v);
}

int ASDPlasticDamageConcrete3DMaterial::setTrialStrainIncr(const Vector& v)
{
	static Vector aux(6);
	// 'strain' is in TENSOR components while the increment arrives in
	// engineering ones, so the sum has to happen in the caller's convention
	aux(0) = strain(0) + v(0);
	aux(1) = strain(1) + v(1);
	aux(2) = strain(2) + v(2);
	aux(3) = 2.0 * strain(3) + v(3);
	aux(4) = 2.0 * strain(4) + v(4);
	aux(5) = 2.0 * strain(5) + v(5);
	return setTrialStrain(aux);
}

int ASDPlasticDamageConcrete3DMaterial::setTrialStrainIncr(const Vector& v, const Vector& /*r*/)
{
	return setTrialStrainIncr(v);
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getStrain(void)
{
	// back to the ENGINEERING convention the caller gave.
	//
	// A MEMBER AND NOT A FUNCTION-LEVEL static, unlike the response getters
	// below. Those are read through the Response machinery, which copies the
	// value out immediately; getStrain() is part of the NDMaterial interface and
	// an element is free to hold the reference while it walks to the next
	// integration point - at which point a shared static would hand it that
	// point's strain instead
	strain_eng(0) = strain(0);
	strain_eng(1) = strain(1);
	strain_eng(2) = strain(2);
	strain_eng(3) = 2.0 * strain(3);
	strain_eng(4) = 2.0 * strain(4);
	strain_eng(5) = 2.0 * strain(5);
	return strain_eng;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getStress(void)
{
	// no conversion: a stress has no engineering convention
	return stress;
}

int ASDPlasticDamageConcrete3DMaterial::commitState(void)
{
	if (implex) {
		// what the extrapolated step delivered is already in stress_implex,
		// recorded by setTrialStrain, and integrate() below does not touch it
		// THE IMPLICIT SOLUTION IS WHAT GETS COMMITTED: the explicit one carried
		// the step, this one carries the state. So this same call measures the
		// error and does the first half of the commit - there is no second
		// return mapping anywhere
		integrate(strain);
		implex_error = implexStressGap(stress_implex, stress);
		// FREEZE what the next explicit pass will use. The direction is the one
		// the return path actually TOOK over the step, not the flow direction at
		// its end: they differ by the curvature of the return path, and the
		// first one is what makes the explicit update reproduce the implicit one
		// EXACTLY when the extrapolated multiplier happens to be right - which
		// is the property the error metric should be measuring the absence of,
		// not something the freezing should introduce on its own
		double dl = dlambda;
		if (dl > 0.0) {
			for (int i = 0; i < 6; ++i) {
				m_commit(i) = (ep_pl(i) - ep_pl_commit(i)) / dl;
				// ep_cr_commit still holds the PREVIOUS committed value here:
				// integrate() reads it and never writes it
				mcr_commit(i) = (ep_cr(i) - ep_cr_commit(i)) / dl;
			}
			kt_rate_commit = (kt - kt_commit) / dl;
			kc_rate_commit = (kc - kc_commit) / dl;
			// the SPLIT ratio is frozen like everything else this pass touches:
			// it is an on/off state of the step
			r_commit = splitWeight();
		}
		else {
			m_commit.Zero();
			mcr_commit.Zero();
			kt_rate_commit = 0.0;
			kc_rate_commit = 0.0;
			r_commit = 0.0;
		}
		dlambda_commit = dl;
		dtime_n_commit = dtime_n;
	}
	strain_commit = strain;
	stress_commit = stress;
	ep_cr_commit = ep_cr;
	ep_pl_commit = ep_pl;
	kt_commit = kt;
	kc_commit = kc;
	// the decomposition of the state just committed, which is what the next
	// explicit pass freezes while letting the reductions move on it
	d_split_commit = d_split;
	V_split_commit = V_split;
	commit_done = true;
	return 0;
}

int ASDPlasticDamageConcrete3DMaterial::revertToLastCommit(void)
{
	strain = strain_commit;
	stress = stress_commit;
	ep_cr = ep_cr_commit;
	ep_pl = ep_pl_commit;
	kt = kt_commit;
	kc = kc_commit;
	d_split = d_split_commit;
	V_split = V_split_commit;
	dtime_n = dtime_n_commit;
	// W IS REBUILT AND NOT STORED. It is a function of the committed split and
	// the two committed measures, through the same expression effective() used
	// to build it at commit time, so the rebuild returns the same operator to
	// the bit and 36 doubles per Gauss point are not carried around for it
	rebuildOperator(d_split_commit, V_split_commit, kt_commit, kc_commit);
	return 0;
}

int ASDPlasticDamageConcrete3DMaterial::revertToStart(void)
{
	strain.Zero();
	strain_commit.Zero();
	stress.Zero();
	stress_commit.Zero();
	stress_implex.Zero();
	sbar.Zero();
	sbar_pos.Zero();
	sbar_neg.Zero();
	ep_cr.Zero();
	ep_cr_commit.Zero();
	ep_pl.Zero();
	ep_pl_commit.Zero();
	kt = kc = 0.0;
	kt_commit = kc_commit = 0.0;
	dlambda = dlambda_commit = 0.0;
	m_commit.Zero();
	mcr_commit.Zero();
	kt_rate_commit = kc_rate_commit = 0.0;
	r_commit = 0.0;
	dtime_n = dtime_n_commit = dtime_0 = 0.0;
	dtime_is_user_defined = false;
	commit_done = false;
	implex_error = 0.0;
	n_iter = 0;
	n_residual = 0;
	bisected = 0;
	plastic = false;
	failed = false;
	stagnated = false;
	extrapolated = false;
	// the split back to EVEN and the operator to the identity: zero eigenvalues
	// put the whole projector in the shared remainder, so PT = PC = I/2
	d_split.Zero();
	d_split_commit.Zero();
	V_split.Zero();
	V_split_commit.Zero();
	for (int i = 0; i < 3; ++i) {
		V_split(i, i) = 1.0;
		V_split_commit(i, i) = 1.0;
	}
	W.Zero();
	for (int i = 0; i < 6; ++i)
		W(i, i) = 1.0;
	C = getInitialTangent();
	return 0;
}

// ==================================================================== //
// 10  the tangent                                                      //
// ==================================================================== //

const Matrix& ASDPlasticDamageConcrete3DMaterial::getInitialTangent(void)
{
	static Matrix D(6, 6);
	double mu2 = E / (1.0 + nu);
	double lam = nu * mu2 / (1.0 - 2.0 * nu);
	elasticMatrix(lam, mu2, D);
	// TO THE ENGINEERING-INPUT CONVENTION: halve the shear COLUMNS. D acts on
	// tensor components internally and OpenSees hands it engineering ones, so
	// this is the second and last place the two conventions meet
	for (int i = 0; i < 6; ++i)
		for (int j = 3; j < 6; ++j)
			D(i, j) *= 0.5;
	return D;
}

const Matrix& ASDPlasticDamageConcrete3DMaterial::getTangent(void)
{
	return C;
}

void ASDPlasticDamageConcrete3DMaterial::computeTangent(void)
{
	// W:C elastic, W:C - A (x) (n:W:C)/den plastic.
	//
	// THE ROTATION OF THE SPLIT PROJECTORS IS NOT DIFFERENTIATED - the same
	// convention the rest of the family uses for a spectral operator, and the
	// reason a numerical tangent disagrees at the states where two eigenvalues
	// meet.
	//
	// UNDER IMPL-EX THIS IS NOT AN APPROXIMATION OF ANYTHING: the explicit
	// update makes sigma affine in eps over the step, so W:C IS the algorithmic
	// tangent - and that is what makes the non-associated flow rule, which
	// normally costs a non-symmetric tangent because m != n is exactly what a
	// user-chosen dilation angle MEANS, cost nothing at the global level.
	static Matrix Cd(6, 6);
	static Matrix WC(6, 6);
	static Vector n(6);
	static Vector m(6);
	static Vector A(6);
	static Vector nWC(6);

	double mu2 = E / (1.0 + nu);
	double lam = nu * mu2 / (1.0 - 2.0 * nu);
	elasticMatrix(lam, mu2, Cd);
	WC.addMatrixProduct(0.0, W, Cd, 1.0);

	// THE LEE-FENVES BRANCH DIFFERENTIATES W, AND THE FARIA ONE DOES NOT, and
	// the asymmetry is not an oversight either way. Faria's W moves because its
	// PROJECTORS rotate, and a rotated projector multiplies a stress that is
	// continuous across the rotation, so the omitted term is small - that is the
	// convention stated above. The Lee-Fenves W is a SCALAR that moves because r
	// moves, and the omitted term is
	//
	//     d sigma / d eps += sbar (x) (dw/dr) (dr/dsbar : C)
	//
	// whose size next to the term that IS kept, w*C, is 1/w. In the fully
	// cracked regime that is sixty-five to one.
	//
	// WHAT IT BROKE, measured on the bench before it was added: the tangent at a
	// MIXED state (r = 0.8123) was off by 40.4% of E against a central
	// difference, and with it 4.9e-10. It is invisible on any converged uniaxial
	// state - dr/dl_i is zero when all three principal stresses share a sign, so
	// a state with sbar = (s, 0, 0) has no term at all - and only the driver's
	// INTERMEDIATE iterates ever see it. That is also why a strain-driven test
	// cannot find it: nothing consults a tangent there.
	//
	// Under IMPL-EX r is frozen by design, so the term is identically zero, and
	// that is precisely what makes sigma affine in eps over the step
	if (damage_combination == DC_LeeFenves && !(implex && extrapolated)) {
		double T = std::abs(d_split(0)) + std::abs(d_split(1))
			+ std::abs(d_split(2));
		if (T > 1.0e-14) {
			double P = 0.0;
			for (int j = 0; j < 3; ++j)
				if (d_split(j) > 0.0) P += d_split(j);
			// dr/dsbar, assembled on the eigen-projectors: dr/dl_j is
			// (H(l_j)*T - P*sign(l_j))/T^2 and dl_j/dsbar is v_j (x) v_j
			static Vector drds(6);
			drds.Zero();
			for (int j = 0; j < 3; ++j) {
				double h = (d_split(j) > 0.0) ? 1.0 : 0.0;
				double sg = (d_split(j) > 0.0) ? 1.0
					: ((d_split(j) < 0.0) ? -1.0 : 0.0);
				double c = (h * T - P * sg) / (T * T);
				double v0 = V_split(0, j);
				double v1 = V_split(1, j);
				double v2 = V_split(2, j);
				drds(0) += c * v0 * v0;
				drds(1) += c * v1 * v1;
				drds(2) += c * v2 * v2;
				drds(3) += c * v0 * v1;
				drds(4) += c * v1 * v2;
				drds(5) += c * v0 * v2;
			}
			double r_now = splitWeight();
			double ot = omega(kt, damage_t, hard_t);
			double oc = omega(kc, damage_c, hard_c);
			double g_t, g_c, g_r;
			leeFenvesDerivatives(ot, oc, r_now, 0.0, 0.0, g_t, g_c, g_r);
			// drds : C, contracted on the FIRST index pair: a symmetric tensor
			// held in Voigt needs its shear entries counted twice, which is the
			// same factor splitWeightRate() carries for the same reason
			static Vector dwde(6);
			for (int j = 0; j < 6; ++j) {
				double s = 0.0;
				for (int i = 0; i < 6; ++i)
					s += ((i < 3) ? 1.0 : 2.0) * drds(i) * Cd(i, j);
				dwde(j) = g_r * s;
			}
			// added BEFORE the halving below, so the engineering convention is
			// applied to it exactly once, like every other column
			for (int i = 0; i < 6; ++i)
				for (int j = 0; j < 6; ++j)
					WC(i, j) += sbar(i) * dwde(j);
		}
	}

	// to the engineering-input convention, once: with the shear columns halved,
	// the row vector n:W:C below comes out with no per-index factor at all
	for (int i = 0; i < 6; ++i)
		for (int j = 3; j < 6; ++j)
			WC(i, j) *= 0.5;

	// THE GUARD IS ON 'DID NOT REACH tol', NOT ON THE VERDICT, and the two
	// stopped being the same thing when the verdict moved into stress units -
	// see residualStress(). It used to read 'failed' alone, so relabelling a
	// step from failed to stagnated silently handed the caller a DIFFERENT
	// operator: the elastoplastic one, built from a state that is not on the
	// surface and whose den therefore means little. Measured on the bench before
	// this line was corrected, that alone moved the delivered stress by 18.3 MPa
	// at damage 0.3/0 - a verdict change is meant to be a report, and it was
	// steering the solve. Reading both flags puts the conservative operator back
	// on every step that did not converge, whatever it ends up being called
	if (implex || !plastic || failed || stagnated) {
		C = WC;
		return;
	}

	double qt, dqt, qc, dqc;
	hard_t.evaluate(kt, qt, dqt);
	hard_c.evaluate(kc, qc, dqc);
	double smax = 0.0;
	dFds(stress, qt, qc, n, smax);
	flowDirection(stress, m);
	double r = splitWeight();
	double h_t, h_c;
	hardeningRates(m, r, h_t, h_c);
	double dF_dkt, dF_dkc;
	dFdKappa(smax, qt, qc, dqt, dqc, dF_dkt, dF_dkc);
	corrector(r, m, h_t, h_c, kt, kc, A);
	double den = contract(n, A) - (dF_dkt * h_t + dF_dkc * h_c);
	// THE SIGN, not just the magnitude. A negative den is the tensile apex,
	// where the correction ADDS stiffness in the flow direction instead of
	// removing it and hands the global Newton an operator that is not even
	// positive definite. W:C is the only defensible answer there, and a step
	// that reached it is a step the return mapping resolved by bracketing rather
	// than by this gradient
	if (!std::isfinite(den) || den <= 1.0e-300) {
		C = WC;
		return;
	}
	for (int j = 0; j < 6; ++j) {
		nWC(j) = n(0) * WC(0, j) + n(1) * WC(1, j) + n(2) * WC(2, j)
			+ 2.0 * (n(3) * WC(3, j) + n(4) * WC(4, j) + n(5) * WC(5, j));
	}
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			C(i, j) = WC(i, j) - A(i) * nWC(j) / den;
}

// ==================================================================== //
// 11  output                                                           //
// ==================================================================== //

void ASDPlasticDamageConcrete3DMaterial::Print(OPS_Stream& s, int flag)
{
	s << "ASDPlasticDamageConcrete3D Material, tag: " << this->getTag() << "\n";
}

int ASDPlasticDamageConcrete3DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// 1000 - elasticity & mass
	if (strcmp(argv[0], "E") == 0) {
		param.setValue(E);
		return param.addObject(1000, this);
	}
	if (strcmp(argv[0], "nu") == 0 || strcmp(argv[0], "v") == 0) {
		param.setValue(nu);
		return param.addObject(1001, this);
	}
	if (strcmp(argv[0], "rho") == 0) {
		param.setValue(rho);
		return param.addObject(1002, this);
	}
	// 2000 - time
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
	// 3000 - the two dials, so that a calibration sweep does not need a rebuild
	if (strcmp(argv[0], "damageT") == 0) {
		param.setValue(damage_t);
		return param.addObject(3000, this);
	}
	if (strcmp(argv[0], "damageC") == 0) {
		param.setValue(damage_c);
		return param.addObject(3001, this);
	}
	if (strcmp(argv[0], "dilatancy") == 0) {
		param.setValue(psi * 180.0 / M_PI);
		return param.addObject(3002, this);
	}
	return -1;
}

int ASDPlasticDamageConcrete3DMaterial::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
	case 1000:
		E = info.theDouble;
		return 0;
	case 1001:
		nu = info.theDouble;
		return 0;
	case 1002:
		rho = info.theDouble;
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
	case 3000:
		damage_t = std::min(1.0, std::max(0.0, info.theDouble));
		return 0;
	case 3001:
		damage_c = std::min(1.0, std::max(0.0, info.theDouble));
		return 0;
	case 3002:
		psi = std::abs(info.theDouble) * M_PI / 180.0;
		return 0;
	default:
		return -1;
	}
}

Vector ASDPlasticDamageConcrete3DMaterial::getHardeningLawVector(
	HardeningLawType ltype, HardeningLawPointComponent c) const
{
	const HardeningLaw& law = (ltype == HardeningLawType::Tension) ? ht : hc;
	Vector out(static_cast<int>(law.points().size()));
	for (std::size_t i = 0; i < law.points().size(); ++i) {
		const ASDHardeningLawPoint& p = law.points()[i];
		switch (c) {
		case HardeningLawPointComponent::TotalStrain:
			out(static_cast<int>(i)) = p.x;
			break;
		case HardeningLawPointComponent::EffectiveStress:
			out(static_cast<int>(i)) = p.q;
			break;
		default:
			out(static_cast<int>(i)) = p.y;
			break;
		}
	}
	return out;
}

Vector ASDPlasticDamageConcrete3DMaterial::getCDPCurveVector(
	HardeningLawType ltype, int component) const
{
	const ASDCDPHardeningCurve& h =
		(ltype == HardeningLawType::Tension) ? hard_t : hard_c;
	Vector out(static_cast<int>(h.size()));
	for (std::size_t i = 0; i < h.size(); ++i) {
		double v = (component == 0) ? h.kappa()[i]
			: ((component == 1) ? h.q()[i] : h.y()[i]);
		out(static_cast<int>(i)) = v;
	}
	return out;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getKappa() const
{
	static Vector d(2);
	d(0) = kt;
	d(1) = kc;
	return d;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getDamage() const
{
	// (1 - omega_t, 1 - omega_c). Identically zero on a side whose dial is 0 -
	// there the whole inelastic strain is plastic and nothing is reduced
	static Vector d(2);
	d(0) = 1.0 - omega(kt, damage_t, hard_t);
	d(1) = 1.0 - omega(kc, damage_c, hard_c);
	return d;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getOmega() const
{
	static Vector d(2);
	d(0) = omega(kt, damage_t, hard_t);
	d(1) = omega(kc, damage_c, hard_c);
	return d;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getSplitWeight() const
{
	// r, Lee and Fenves' tensile weight of the CURRENT effective stress.
	//
	// Published because it was the one quantity of the formulation that nothing
	// could see from outside, and it decides three separate things: how the flow
	// is split between the permanent strain and the reductions, how the two
	// hardening measures advance, and - under '-damageCombination leeFenves' -
	// how the two reductions combine. It is also the quantity the return
	// mapping's limit cycle swings on, so a trace that cannot read it cannot
	// diagnose that failure.
	//
	// In the Lee-Fenves branch this makes the scalar reduction fully observable
	// from published quantities: with 'omega' and the two weights it reproduces
	// 1-d exactly, and it also has to equal the ratio of 'stress' to
	// 'effectiveStress', which is an independent check of the same number
	static Vector d(1);
	d(0) = splitWeight();
	return d;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getStrength() const
{
	static Vector d(2);
	double dq;
	hard_t.evaluate(kt, d(0), dq);
	hard_c.evaluate(kc, d(1), dq);
	return d;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getPlasticStrainVector() const
{
	// back to the ENGINEERING convention, as a strain-like output must be
	static Vector d(6);
	d(0) = ep_pl(0);
	d(1) = ep_pl(1);
	d(2) = ep_pl(2);
	d(3) = 2.0 * ep_pl(3);
	d(4) = 2.0 * ep_pl(4);
	d(5) = 2.0 * ep_pl(5);
	return d;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getCrackingStrainVector() const
{
	// the share of the flow the REDUCTIONS carry rather than the plastic strain.
	// With the plastic strain it sums to lambda*m exactly, which is the
	// bookkeeping invariant
	static Vector d(6);
	d(0) = ep_cr(0);
	d(1) = ep_cr(1);
	d(2) = ep_cr(2);
	d(3) = 2.0 * ep_cr(3);
	d(4) = 2.0 * ep_cr(4);
	d(5) = 2.0 * ep_cr(5);
	return d;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getEffectiveStress() const
{
	return sbar;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getImplexError() const
{
	static Vector d(1);
	d(0) = implex_error;
	return d;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getImplexStress() const
{
	return stress_implex;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getTimeIncrements() const
{
	static Vector d(3);
	d(0) = dtime_n;
	d(1) = dtime_n_commit;
	d(2) = dtime_0;
	return d;
}

const Vector& ASDPlasticDamageConcrete3DMaterial::getIntegrationInfo() const
{
	static Vector d(7);
	d(0) = static_cast<double>(n_iter);
	d(1) = static_cast<double>(n_residual);
	d(2) = dlambda;
	d(3) = plastic ? 1.0 : 0.0;
	d(4) = failed ? 1.0 : 0.0;
	d(5) = stagnated ? 1.0 : 0.0;
	d(6) = static_cast<double>(bisected);
	return d;
}

Response* ASDPlasticDamageConcrete3DMaterial::setResponse(
	const char** argv, int argc, OPS_Stream& output)
{
	auto make_resp = [&output, this](int rid, const Vector& v,
		const std::vector<std::string>* labels = nullptr) -> MaterialResponse* {
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

	static std::vector<std::string> lb_pair_tc = { "T", "C" };
	static std::vector<std::string> lb_damage = { "d+", "d-" };
	static std::vector<std::string> lb_omega = { "w+", "w-" };
	static std::vector<std::string> lb_split_weight = { "r" };
	static std::vector<std::string> lb_tensor = { "11", "22", "33", "12", "23", "13" };
	static std::vector<std::string> lb_implex_error = { "Error" };
	static std::vector<std::string> lb_time = { "dTime", "dTimeCommit", "dTimeInitial" };
	static std::vector<std::string> lb_info = { "nIter", "nResidual", "lambda",
		"plastic", "failed", "stagnated", "bisected" };

	if (argc > 0) {
		// 1000 - the INPUT backbones, as given (and regularized)
		if (strcmp(argv[0], "Ce") == 0)
			return make_resp(1000, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::TotalStrain));
		if (strcmp(argv[0], "Cs") == 0)
			return make_resp(1001, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::NominalStress));
		if (strcmp(argv[0], "Te") == 0)
			return make_resp(1002, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::TotalStrain));
		if (strcmp(argv[0], "Ts") == 0)
			return make_resp(1003, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::NominalStress));
		// 1100 - the CDP hardening curves the surface actually sees. Worth
		// exposing separately from the input: the conversion kappa = x - q/E is
		// where a table stops being interpolable, and a user looking at a
		// surprising response needs to see the curve the surface saw.
		//
		// THE PROTOTYPE IS SIZED OFF THE LAW AND NOT OFF THE CURVE, because the
		// curve does not exist yet: its conversion is lazy - it needs the
		// element's characteristic length - and a recorder is set up before the
		// first setTrialStrain. Calling prepare() here would be worse than an
		// empty prototype: it would freeze lch at 1 before the element could
		// supply it, and silently disable the regularization. The law's point
		// count is an upper bound on the curve's, and Vector::operator= resizes
		// on the first read anyway
		if (strcmp(argv[0], "Tkappa") == 0)
			return make_resp(1100, Vector(static_cast<int>(ht.points().size())));
		if (strcmp(argv[0], "Tq") == 0)
			return make_resp(1101, Vector(static_cast<int>(ht.points().size())));
		if (strcmp(argv[0], "Ckappa") == 0)
			return make_resp(1102, Vector(static_cast<int>(hc.points().size())));
		if (strcmp(argv[0], "Cq") == 0)
			return make_resp(1103, Vector(static_cast<int>(hc.points().size())));
		// 2000 - state
		if (strcmp(argv[0], "kappa") == 0 || strcmp(argv[0], "Kappa") == 0 ||
			strcmp(argv[0], "equivalentPlasticStrain") == 0 ||
			strcmp(argv[0], "EquivalentPlasticStrain") == 0)
			return make_resp(2000, getKappa(), &lb_pair_tc);
		if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0)
			return make_resp(2001, getDamage(), &lb_damage);
		if (strcmp(argv[0], "omega") == 0 || strcmp(argv[0], "Omega") == 0)
			return make_resp(2002, getOmega(), &lb_omega);
		if (strcmp(argv[0], "strength") == 0 || strcmp(argv[0], "Strength") == 0)
			return make_resp(2003, getStrength(), &lb_pair_tc);
		if (strcmp(argv[0], "splitWeight") == 0 || strcmp(argv[0], "SplitWeight") == 0)
			return make_resp(2007, getSplitWeight(), &lb_split_weight);
		if (strcmp(argv[0], "plasticStrain") == 0 || strcmp(argv[0], "PlasticStrain") == 0)
			return make_resp(2004, getPlasticStrainVector(), &lb_tensor);
		if (strcmp(argv[0], "crackingStrain") == 0 || strcmp(argv[0], "CrackingStrain") == 0)
			return make_resp(2005, getCrackingStrainVector(), &lb_tensor);
		if (strcmp(argv[0], "effectiveStress") == 0 || strcmp(argv[0], "EffectiveStress") == 0)
			return make_resp(2006, getEffectiveStress(), &lb_tensor);
		// 3000 - IMPL-EX and the integration
		if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0)
			return make_resp(3000, getImplexError(), &lb_implex_error);
		// the stress the step DELIVERED, which after the commit is otherwise
		// gone: under IMPL-EX the implicit re-solve installs its own over it
		if (strcmp(argv[0], "implexStress") == 0 || strcmp(argv[0], "ImplexStress") == 0)
			return make_resp(3003, getImplexStress(), &lb_tensor);
		if (strcmp(argv[0], "timeIncrements") == 0 || strcmp(argv[0], "TimeIncrements") == 0)
			return make_resp(3001, getTimeIncrements(), &lb_time);
		if (strcmp(argv[0], "integrationInfo") == 0 || strcmp(argv[0], "IntegrationInfo") == 0)
			return make_resp(3002, getIntegrationInfo(), &lb_info);
	}
	return NDMaterial::setResponse(argv, argc, output);
}

int ASDPlasticDamageConcrete3DMaterial::getResponse(int responseID, Information& matInformation)
{
	switch (responseID) {
	case 1000: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::TotalStrain));
	case 1001: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::NominalStress));
	case 1002: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::TotalStrain));
	case 1003: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::NominalStress));
	case 1100: return matInformation.setVector(getCDPCurveVector(HardeningLawType::Tension, 0));
	case 1101: return matInformation.setVector(getCDPCurveVector(HardeningLawType::Tension, 1));
	case 1102: return matInformation.setVector(getCDPCurveVector(HardeningLawType::Compression, 0));
	case 1103: return matInformation.setVector(getCDPCurveVector(HardeningLawType::Compression, 1));
	case 2000: return matInformation.setVector(getKappa());
	case 2001: return matInformation.setVector(getDamage());
	case 2002: return matInformation.setVector(getOmega());
	case 2003: return matInformation.setVector(getStrength());
	case 2004: return matInformation.setVector(getPlasticStrainVector());
	case 2005: return matInformation.setVector(getCrackingStrainVector());
	case 2006: return matInformation.setVector(getEffectiveStress());
	case 2007: return matInformation.setVector(getSplitWeight());
	case 3000: return matInformation.setVector(getImplexError());
	case 3003: return matInformation.setVector(getImplexStress());
	case 3001: return matInformation.setVector(getTimeIncrements());
	case 3002: return matInformation.setVector(getIntegrationInfo());
	default:
		break;
	}
	return NDMaterial::getResponse(responseID, matInformation);
}

// ==================================================================== //
// 12  serialization                                                    //
// ==================================================================== //

int ASDPlasticDamageConcrete3DMaterial::sendSelf(int commitTag, Channel& theChannel)
{
	int counter;

	// THE CDP CURVES ARE NOT SENT: they are a pure function of the laws, of E and
	// of the damage cap, so the receiver rebuilds them. What has to travel is
	// regularization_done, or the receiver would regularize laws that already are
	int nv_dbl = 180 + ht.serializationDataSize() + hc.serializationDataSize();

	static ID idata(13);
	counter = 0;
	idata(counter++) = getTag();
	idata(counter++) = static_cast<int>(implex);
	idata(counter++) = static_cast<int>(implex_control);
	idata(counter++) = static_cast<int>(implex_abort_on_error);
	idata(counter++) = static_cast<int>(auto_regularize);
	idata(counter++) = static_cast<int>(regularization_done);
	idata(counter++) = static_cast<int>(dtime_is_user_defined);
	idata(counter++) = static_cast<int>(commit_done);
	idata(counter++) = max_iter;
	idata(counter++) = max_backtrack;
	idata(counter++) = static_cast<int>(damage_combination);
	idata(counter++) = (static_cast<int>(plastic))
		| (static_cast<int>(failed) << 1)
		| (static_cast<int>(stagnated) << 2)
		| (static_cast<int>(extrapolated) << 3);
	idata(counter++) = nv_dbl;
	if (theChannel.sendID(getDbTag(), commitTag, idata) < 0) {
		opserr << "ASDPlasticDamageConcrete3DMaterial::sendSelf() - failed to send INT data\n";
		return -1;
	}

	Vector ddata(nv_dbl);
	counter = 0;
	ddata(counter++) = E;
	ddata(counter++) = nu;
	ddata(counter++) = rho;
	ddata(counter++) = psi;
	ddata(counter++) = ecc;
	ddata(counter++) = fb0_fc0;
	ddata(counter++) = Kc;
	ddata(counter++) = alpha;
	ddata(counter++) = gam;
	ddata(counter++) = damage_t;
	ddata(counter++) = damage_c;
	ddata(counter++) = stiffness_recovery_t;
	ddata(counter++) = stiffness_recovery_c;
	ddata(counter++) = implex_error_tolerance;
	ddata(counter++) = implex_time_redution_limit;
	ddata(counter++) = implex_alpha;
	ddata(counter++) = tol;
	ddata(counter++) = stagnation_tol;
	ddata(counter++) = lch;
	ddata(counter++) = lch_ref;
	ddata(counter++) = kt;
	ddata(counter++) = kc;
	ddata(counter++) = kt_commit;
	ddata(counter++) = kc_commit;
	ddata(counter++) = dlambda;
	ddata(counter++) = dlambda_commit;
	ddata(counter++) = kt_rate_commit;
	ddata(counter++) = kc_rate_commit;
	ddata(counter++) = r_commit;
	ddata(counter++) = dtime_n;
	ddata(counter++) = dtime_n_commit;
	ddata(counter++) = dtime_0;
	ddata(counter++) = implex_error;
	ddata(counter++) = static_cast<double>(n_iter);
	ddata(counter++) = static_cast<double>(n_residual);
	ddata(counter++) = static_cast<double>(bisected);
	for (int i = 0; i < 6; ++i) ddata(counter++) = strain(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = strain_commit(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = stress(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = stress_commit(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = stress_implex(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = sbar(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = sbar_pos(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = sbar_neg(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = ep_cr(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = ep_cr_commit(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = ep_pl(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = ep_pl_commit(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = m_commit(i);
	for (int i = 0; i < 6; ++i) ddata(counter++) = mcr_commit(i);
	for (int i = 0; i < 3; ++i) ddata(counter++) = d_split(i);
	for (int i = 0; i < 3; ++i) ddata(counter++) = d_split_commit(i);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			ddata(counter++) = V_split(i, j);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			ddata(counter++) = V_split_commit(i, j);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			ddata(counter++) = W(i, j);
	ht.serialize(ddata, counter);
	hc.serialize(ddata, counter);
	if (counter != nv_dbl) {
		opserr << "ASDPlasticDamageConcrete3DMaterial::sendSelf() - internal error: "
			"wrote " << counter << " doubles, declared " << nv_dbl << "\n";
		return -1;
	}
	if (theChannel.sendVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDPlasticDamageConcrete3DMaterial::sendSelf() - failed to send DBL data\n";
		return -1;
	}
	return 0;
}

int ASDPlasticDamageConcrete3DMaterial::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int counter;

	static ID idata(13);
	if (theChannel.recvID(getDbTag(), commitTag, idata) < 0) {
		opserr << "ASDPlasticDamageConcrete3DMaterial::recvSelf() - failed to receive INT data\n";
		return -1;
	}
	counter = 0;
	setTag(idata(counter++));
	implex = static_cast<bool>(idata(counter++));
	implex_control = static_cast<bool>(idata(counter++));
	implex_abort_on_error = static_cast<bool>(idata(counter++));
	auto_regularize = static_cast<bool>(idata(counter++));
	regularization_done = static_cast<bool>(idata(counter++));
	dtime_is_user_defined = static_cast<bool>(idata(counter++));
	commit_done = static_cast<bool>(idata(counter++));
	max_iter = idata(counter++);
	max_backtrack = idata(counter++);
	damage_combination = static_cast<DamageCombination>(idata(counter++));
	int flags = idata(counter++);
	plastic = (flags & 1) != 0;
	failed = (flags & 2) != 0;
	stagnated = (flags & 4) != 0;
	extrapolated = (flags & 8) != 0;
	int nv_dbl = idata(counter++);

	Vector ddata(nv_dbl);
	if (theChannel.recvVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDPlasticDamageConcrete3DMaterial::recvSelf() - failed to receive DBL data\n";
		return -1;
	}
	counter = 0;
	E = ddata(counter++);
	nu = ddata(counter++);
	rho = ddata(counter++);
	psi = ddata(counter++);
	ecc = ddata(counter++);
	fb0_fc0 = ddata(counter++);
	Kc = ddata(counter++);
	alpha = ddata(counter++);
	gam = ddata(counter++);
	damage_t = ddata(counter++);
	damage_c = ddata(counter++);
	stiffness_recovery_t = ddata(counter++);
	stiffness_recovery_c = ddata(counter++);
	implex_error_tolerance = ddata(counter++);
	implex_time_redution_limit = ddata(counter++);
	implex_alpha = ddata(counter++);
	tol = ddata(counter++);
	stagnation_tol = ddata(counter++);
	lch = ddata(counter++);
	lch_ref = ddata(counter++);
	kt = ddata(counter++);
	kc = ddata(counter++);
	kt_commit = ddata(counter++);
	kc_commit = ddata(counter++);
	dlambda = ddata(counter++);
	dlambda_commit = ddata(counter++);
	kt_rate_commit = ddata(counter++);
	kc_rate_commit = ddata(counter++);
	r_commit = ddata(counter++);
	dtime_n = ddata(counter++);
	dtime_n_commit = ddata(counter++);
	dtime_0 = ddata(counter++);
	implex_error = ddata(counter++);
	n_iter = static_cast<int>(ddata(counter++));
	n_residual = static_cast<int>(ddata(counter++));
	bisected = static_cast<int>(ddata(counter++));
	for (int i = 0; i < 6; ++i) strain(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) strain_commit(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) stress(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) stress_commit(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) stress_implex(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) sbar(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) sbar_pos(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) sbar_neg(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) ep_cr(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) ep_cr_commit(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) ep_pl(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) ep_pl_commit(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) m_commit(i) = ddata(counter++);
	for (int i = 0; i < 6; ++i) mcr_commit(i) = ddata(counter++);
	for (int i = 0; i < 3; ++i) d_split(i) = ddata(counter++);
	for (int i = 0; i < 3; ++i) d_split_commit(i) = ddata(counter++);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			V_split(i, j) = ddata(counter++);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			V_split_commit(i, j) = ddata(counter++);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			W(i, j) = ddata(counter++);
	ht.deserialize(ddata, counter);
	hc.deserialize(ddata, counter);

	// rebuild what was not sent. curves_ready is deliberately left false so that
	// prepare() converts the (already regularized) laws on the first
	// setTrialStrain, and the tangent is put back in a usable state right away
	curves_ready = false;
	setup_reported = false;
	prepare();
	computeTangent();
	return 0;
}
