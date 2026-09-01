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
// A J2 (von Mises) plasticity model for steel with non-linear kinematic
// hardening (Armstrong-Frederick / Chaboche, two backstresses) and an IMPL-EX
// integration scheme.
//
// THE 3D SIBLING OF ASDSteel1DMaterial, and deliberately only two thirds of it:
// the elasto-plastic core and IMPL-EX. The 1D's fracture (a scalar damage keyed
// to an empirically regularized ultimate strain) and its buckling RVE are NOT
// here - the second one because a buckling bar is a one-dimensional object and
// the notion does not survive the move to a continuum point.
//
// WHAT MAKES IT THE SAME MODEL. With
//
//     f     = ||s - a|| - sqrt(2/3)*sy          a = a1 + a2, sy CONSTANT
//     dep   = dgamma * n                        n = (s - a)/||s - a||
//     dp    = sqrt(2/3) * dgamma
//     da_i  = (2/3)*H_i*dep - gamma_i*a_i*dp
//
// the uniaxial-stress restriction is EXACTLY the 1D model, with the same H_i and
// the same gamma_i - no rescaling. The proof is short enough to keep here.
// Uniaxially s = sig*diag(2/3,-1/3,-1/3) and a_i = A_i*diag(1,-1/2,-1/2), and
// since ||diag(1,-1/2,-1/2)|| = sqrt(3/2),
//
//     ||s - a|| = sqrt(2/3) * |sig - (3/2)*A|      so  alpha_i = (3/2)*A_i
//     n11       = sg*sqrt(2/3)                     so  dep11 = sg*dp
//     dA_i      = (2/3)*H_i*sg*dp - gamma_i*A_i*dp
//
// and multiplying the last line by 3/2 gives d(alpha_i) = (H_i*sg -
// gamma_i*alpha_i)*dp, which is the 1D law with dp playing the part of the 1D
// lambda. Two consequences worth stating: the saturated ||a_i|| is
// sqrt(2/3)*H_i/gamma_i, so the uniaxial asymptote is sy + sum(H_i/gamma_i) =
// su, which is why stressReference() below is IDENTICAL to the 1D one and an
// IMPL-EX tolerance means the same thing in both; and pure shear yields at
// sy/sqrt(3), which is the cheapest independent check of every factor of two in
// the file.
//
// THE RETURN MAPPING FREEZES THE FLOW DIRECTION at the trial relative stress,
// which turns the whole thing into a scalar Newton on one multiplier and keeps
// the backstresses on the same exact exponential map the 1D uses. It is exact
// for any radial path in relative-stress space - all uniaxial loading, monotonic
// or cyclic, and pure shear - and first order otherwise, which is the same
// formal order a fully implicit closest-point projection has. What freezing buys
// is one unknown instead of seven; what it costs is the error constant on
// strongly non-proportional paths.
//
// VOIGT CONVENTION, and it is the same rule ASDPlasticDamageConcrete3DMaterial
// states: EVERY 6-VECTOR INSIDE THIS FILE CARRIES THE TENSOR COMPONENTS OF THE
// SHEAR. The engineering convention that OpenSees speaks lives in
// setTrialStrain, in getStrain, in the two tangents and in the plastic-strain
// response, and nowhere else. Inside, that makes the plastic update a plain
// 'ep += dgamma*n' with no weight matrix anywhere, and leaves the factor two in
// exactly one place: the full contraction of two stress-like tensors.

#ifndef ASDSteel3DMaterial_h
#define ASDSteel3DMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <IMPLEXManager.h>
#include <cmath>

// This material runs an IMPL-EX scheme, so it takes part in the IMPL-EX error
// control: it registers itself, says when it took part in a step, and measures
// its own extrapolation error when asked. See IMPLEXManager.h
class ASDSteel3DMaterial : public NDMaterial, public IMPLEXObject
{
public:
	/**
	How the tangent is built on the IMPLICIT path. Under IMPL-EX none of this
	applies: there the elastic operator IS the algorithmic tangent, exactly and
	not approximately, because the extrapolated update is affine in the strain.
	*/
	enum TangentType {
		// THE DEFAULT. The consistent tangent with the AF memory term taken
		// collinear with the frozen normal. Symmetric, and EXACT on any radial
		// path - which is every verification case and most structural ones.
		//
		// The default is not about accuracy, it is about the solver: the full
		// tangent is non-symmetric, and a symmetric OpenSees solver symmetrizes
		// whatever it is handed, in silence. A tangent that cannot be delivered
		// is not an improvement.
		Tangent_Collinear = 0,
		// The full consistent tangent, including the non-collinearity of the AF
		// memory term. NON-SYMMETRIC - and not because of the frozen normal:
		// the Armstrong-Frederick recovery term is non-associative, so every
		// consistent tangent of this model is. Needs a non-symmetric solver to
		// be worth anything.
		Tangent_Full = 1,
		// Six perturbed solves. For verification, and as the way out if the
		// analytical one is ever found wrong. Ignored under IMPL-EX.
		Tangent_Numerical = 2,
		// The elastic operator. Always admissible, never quadratic.
		Tangent_Elastic = 3
	};

	/**
	Everything a user can type, plus what the parser derives from it.
	*/
	class InputParameters {
	public:
		// elasticity and mass
		double E = 0.0;
		double nu = 0.0;
		double rho = 0.0;
		// yield stress. The VIRGIN radius: with the Bauschinger dial at zero it
		// is constant, exactly as in the 1D; with the dial on, the identity
		// pairs below contract it with the accumulated plastic strain
		double sy = 0.0;
		// Chaboche / Armstrong-Frederick kinematic hardening, two backstresses.
		// Calibrated by the parser from (E, sy, su, eu) with the same formula
		// the 1D uses, so that sy + H1/gamma1 + H2/gamma2 == su exactly
		double H1 = 0.0;
		double gamma1 = 0.0;
		double H2 = 0.0;
		double gamma2 = 0.0;
		// THE BAUSCHINGER DIAL, in [0, 1]: 0 keeps the original calibration
		// (bit-identical response), 1 is the recipe tuned on Menegotto-Pinto
		// (OpenSees-Testing/asd-steel-3d, ex04/ex05). Each identity pair is a
		// Voce softening term (amplitude pair_Q >= 0, rate pair_b) that ships
		// with a kinematic twin at the SAME rate (H = Q*b, gamma = b), derived
		// - never stored - so an unpaired softening cannot exist. On the
		// monotonic backbone the twin cancels its softening identically, so
		// the backbone is the original one for EVERY dial value; on reversal
		// the twin must travel 2*Q to change sign while the softening only
		// follows the accumulated p - the released difference IS the
		// Bauschinger effect. The twin also protects the return mapping: it
		// steepens the scalar Newton residual two to one where the softening
		// flattens it, so the monotonicity margin never drops below the
		// elastic unit. Same numbers, same names, same meaning as in
		// ASDSteel1DMaterial::InputParameters - the isotropic term is
		// IDENTICAL in the two models (f = ||xi|| - sqrt(2/3)*(sy + R(p)), and
		// the 3D p maps one-to-one onto the 1D lambda), so the uniaxial
		// equivalence holds with no new proof
		double bauschinger = 0.0;
		static constexpr int NPAIRS = 3;
		static constexpr int NKIN = 2 + NPAIRS;
		double pair_Q[NPAIRS] = { 0.0, 0.0, 0.0 };
		double pair_b[NPAIRS] = { 0.0, 0.0, 0.0 };
		inline bool hasBauschinger() const { return bauschinger > 0.0; }
		// the isotropic part of the yield radius (uniaxial measure, the
		// sqrt(2/3) is applied where sy gets its own), and its derivative
		inline double isoR(double p) const {
			double r = 0.0;
			for (int j = 0; j < NPAIRS; ++j)
				r -= pair_Q[j] * (1.0 - std::exp(-pair_b[j] * p));
			return r;
		}
		inline double isoRprime(double p) const {
			double r = 0.0;
			for (int j = 0; j < NPAIRS; ++j)
				r -= pair_Q[j] * pair_b[j] * std::exp(-pair_b[j] * p);
			return r;
		}
		// the kinematic terms, originals first, twins after. Fills H and g
		// (sized NKIN) and returns how many are active: 2 when the dial is at
		// zero, so the default pays nothing for the machinery
		inline int kinTerms(double* H, double* g) const {
			H[0] = H1; g[0] = gamma1;
			H[1] = H2; g[1] = gamma2;
			if (!hasBauschinger())
				return 2;
			for (int j = 0; j < NPAIRS; ++j) {
				H[2 + j] = pair_Q[j] * pair_b[j];
				g[2 + j] = pair_b[j];
			}
			return NKIN;
		}
		// IMPL-EX
		bool implex = false;
		bool implex_control = false;
		// LEGACY: let the material itself fail the step when its own error is
		// out of tolerance. Off by default, and it should stay off: a material
		// that aborts owns a policy it cannot see the analysis to choose. The
		// rejection belongs to the convergence test wrapper - see IMPLEXManager.h
		bool implex_abort_on_error = false;
		double implex_error_tolerance = 0.05;
		double implex_time_redution_limit = 0.01;
		double implex_alpha = 1.0;
		// how the implicit tangent is built
		int tangent_type = Tangent_Collinear;
	};

	/**
	Everything the implicit pass of integrate() can write. Saved and restored
	around the non-destructive measurement of the IMPL-EX error, so that a step
	that is measured and then rejected is a step that never happened.

	It is SMALL, and deliberately so. The frozen flow direction n_commit - the
	3D counterpart of the 1D's sg_commit, and the one committed quantity an
	implicit pass has any business overwriting - is written in commitState()
	and NOT in integrate(), so it never needs putting back. That is the whole
	reason the freezing lives where it lives: in the 1D it happens inside
	compute(), which is why that material needs sg_commit in its snapshot.

	C is in here because the recorders run between convergence and commit and
	must see the operator the step delivered, not the one the measurement built.
	*/
	struct TrialState {
		Vector stress = Vector(6);
		Vector ep = Vector(6);
		Vector a[InputParameters::NKIN] = { Vector(6), Vector(6), Vector(6), Vector(6), Vector(6) };
		Vector n_flow = Vector(6);
		Matrix C = Matrix(6, 6);
		double p = 0.0;
		double dgamma = 0.0;
		bool plastic = false;
		bool failed = false;
	};

public:
	// life-cycle
	ASDSteel3DMaterial(int _tag, const InputParameters& _params);
	ASDSteel3DMaterial();
	~ASDSteel3DMaterial();

	// info
	const char* getClassType(void) const { return "ASDSteel3DMaterial"; }

	// density
	double getRho(void);

	// set strain
	int setTrialStrain(const Vector& v);
	int setTrialStrain(const Vector& v, const Vector& r);
	int setTrialStrainIncr(const Vector& v);
	int setTrialStrainIncr(const Vector& v, const Vector& r);

	// get state
	const Vector& getStrain(void);
	const Vector& getStress(void);
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);

	// handle state
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// copy and others...
	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* code);
	const char* getType(void) const;
	int getOrder(void) const;
	void Print(OPS_Stream& s, int flag = 0);

	// send/recv self
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);

	// IMPL-EX error control (see IMPLEXManager.h)
	double computeImplexErrorMetric(void);
	double implexTimeRatio(void) const;

private:
	// the implicit answer at the current trial strain, from the committed state.
	// Returns 0, or -1 when the scalar Newton did not converge
	int integrate(void);
	// the extrapolated answer, from the committed state and the frozen direction
	void extrapolate(void);
	// dtime_n / dtime_n_commit * alpha, clamped at zero
	double timeFactor(void) const;
	// the implicit tangent, per params.tangent_type. Under IMPL-EX it is C_e
	void computeTangent(void);
	// the -tangent path: seven solves, and it builds C itself
	int numericalTangent(void);
	// the elastic operator, tensor-in / tensor-out (no engineering halving)
	void elasticTangent(Matrix& out) const;
	// the largest stress this material can carry, the denominator of the metric
	double stressReference(void) const;
	// THE metric: the largest Voigt component of the gap, over stressReference()
	double implexStressGap(const Vector& delivered, const Vector& stress_implicit) const;
	// save/restore everything the implicit pass writes
	void saveTrialState(TrialState& x) const;
	void restoreTrialState(const TrialState& x);
	// responses
	const Vector& getEquivalentPlasticStrain() const;
	const Vector& getBackStress1() const;
	const Vector& getBackStress2() const;
	const Vector& getPlasticStrain() const;
	const Vector& getYieldFunction() const;
	const Vector& getImplexError() const;
	const Vector& getImplexStress() const;
	const Vector& getTimeIncrements() const;

private:
	// input
	InputParameters params;

	// state variables - strain and stress. TENSOR components of the shear, see
	// the note at the top of this file
	Vector strain = Vector(6);
	Vector strain_commit = Vector(6);
	Vector stress = Vector(6);
	Vector stress_commit = Vector(6);
	// WHAT THE STEP DELIVERED TO THE ELEMENT, kept because it is otherwise
	// unobservable: under IMPL-EX commitState re-solves implicitly and installs
	// that answer over 'stress', so from outside - a recorder, a test - the
	// extrapolated stress is gone by the time anyone can ask. Written in ONE
	// place, at the end of setTrialStrain, covering every branch
	Vector stress_implex = Vector(6);

	// state variables - plasticity. The backstresses: the 2 originals first,
	// then the twins of the identity pairs (zero, and untouched, at dial zero)
	Vector ep = Vector(6);
	Vector ep_commit = Vector(6);
	Vector a[InputParameters::NKIN] = { Vector(6), Vector(6), Vector(6), Vector(6), Vector(6) };
	Vector a_commit[InputParameters::NKIN] = { Vector(6), Vector(6), Vector(6), Vector(6), Vector(6) };
	double p = 0.0;
	double p_commit = 0.0;
	double dgamma = 0.0;

	// the flow direction of the current pass, and THE FROZEN ONE the next
	// explicit pass will use. n_commit is zero exactly when the committed step
	// was elastic, which is the 3D counterpart of the 1D's sg_commit == 0 and
	// degrades the extrapolation to a purely elastic prediction
	Vector n_flow = Vector(6);
	Vector n_commit = Vector(6);
	double dgamma_commit = 0.0;

	// tangent
	Matrix C = Matrix(6, 6);

	// state variables - implex
	double dtime_n = 0.0;
	double dtime_n_commit = 0.0;
	double dtime_0 = 0.0;
	bool dtime_is_user_defined = false;
	bool commit_done = false;
	double implex_error = 0.0;

	// diagnostics
	bool plastic = false;
	bool failed = false;

	// what computeTangent() needs from the last implicit pass. Valid only right
	// after a plastic integrate(), never read otherwise, never serialized, and
	// deliberately NOT in TrialState: nothing reads them after the tangent is
	// built, and the peek puts C itself back
	double tg_r = 0.0;      // ||xi_trial||
	double tg_beta = 0.0;   // the component of xi along the frozen normal
	double tg_Rp = 0.0;     // dR/dDp at the solution
	Vector tg_xih = Vector(6);

	// RESPONSE BUFFERS, one per getter and held as members on purpose. Every
	// getter in ASDSteel1D returns a reference to a function-local static, i.e.
	// ONE buffer shared by every material point in the process; it survives only
	// because callers happen to read one at a time. Thirty doubles per gauss
	// point is not a price worth that
	mutable Vector out_ple = Vector(1);
	mutable Vector out_a1 = Vector(6);
	mutable Vector out_a2 = Vector(6);
	mutable Vector out_ep = Vector(6);
	mutable Vector out_yield = Vector(1);
	mutable Vector out_error = Vector(1);
	mutable Vector out_impstress = Vector(6);
	mutable Vector out_time = Vector(3);
};

#endif // ASDSteel3DMaterial_h
