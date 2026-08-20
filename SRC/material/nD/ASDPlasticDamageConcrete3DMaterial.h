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
// Concrete Damaged Plasticity with the DAMAGE IN THE PREDICTION and the
// plasticity on the NOMINAL surface. One yield surface, one plastic multiplier,
// a prescribed dilation angle, and the return mapping run on the stress that
// equilibrium actually sees.
//
// THE ORDER IS THE WHOLE IDEA. The classical CDP puts the plasticity in the
// effective stress and the damage on top of it; here the damage comes FIRST, in
// the elastic prediction, and the return mapping runs on what is left:
//
//     sbar  = C : (eps - eps_p)                 eps_p = the PLASTIC strain
//     sbar  = PT:sbar + PC:sbar                 spectral split, by sign
//     sigma = omega_t*PT:sbar + omega_c*PC:sbar the DAMAGED prediction
//     F(sigma, qt(kt), qc(kc)) <= 0             one surface, NOMINAL backbones
//
// Two consequences, and they are the reason for the change of formulation:
//
// * the hardening curves are the NOMINAL backbones. The classical CDP asks for
//   q = y/(1-d), which GROWS where y softens - 3.0 -> 1.13 -> 9.45 -> 434.85 MPa
//   on the Model Code tensile preset - and ends with points at constant kappa,
//   i.e. an infinite hardening modulus, which is not a curve a return mapping
//   can be handed. Here the surface sees the real stress, so its input is the
//   real curve, and the input damage arrays are not needed at all;
// * sigma : deps_p is the actual dissipation, and the flow direction is measured
//   against the stress the material carries rather than against a fictitious one
//   that keeps growing while the material softens.
//
// THE REDUCTION HAS A DIRECTION, and that is not an option. The third line
// above applies the two reductions to the two SPECTRAL PARTS separately, after
// Faria, so a cracked side loses only its positive part and the compressive one
// is untouched: the recovery of stiffness on crack closure is TOTAL and
// AUTOMATIC, performed by the split itself, with nothing to calibrate.
//
// THE SCALAR ALTERNATIVE WAS BUILT, MEASURED AND REMOVED. It was the classical
// Lee-Fenves combination, sigma = (1-d)*sbar with (1-d) = (1-s_t*d_c)(1-s_c*d_t)
// and the split weight saturated so that a nearly tensile state counts as
// tensile - cheaper per step and, on a monotone path, indistinguishable. What
// decided against it is not its response but its CONTROLLABILITY: under IMPL-EX
// it converges at the same first order with 3.3 times the constant and enters
// the asymptotic regime two levels later, and with the step under error control
// it is the one branch where tightening the tolerance from 1e-2 to 1e-3 made the
// trajectory WORSE - 1.372 to 2.470 MPa at the reversals - by leaving onto a
// more dilated path that the gate cannot see and cannot undo. The measurements
// are in the OpenSees-Testing note asd-new-mats-implex-control/doc/asd_cdp3d.tex
// (ex53/ex54 and their controlled versions); the branch itself is in the history
// of this file.
//
// WHERE EACH UNLOADING LANDS IS A USER PARAMETER, one per side: a DAMAGE FACTOR
// df in [0, 1]. 0 = pure plasticity (elastic unloading onto the full inelastic
// strain), 1 = pure damage (secant unloading onto the origin), and every value
// between reproduces the SAME backbone with residual (1-df)*kappa and secant
// omega*E. The pair (df_t, df_c) = (1, 0) - cracking in tension, crushing in
// compression - is the default. The reduction is DERIVED and not given:
//
//     omega(kappa) = q / (df*E*kappa + q)
//
// which needs no damage array, is 1 at kappa = 0 at any df, and is identically 1
// at df = 0, so the pure-plasticity limit costs not one operation.
//
// THE FLOW IS SHARED between the permanent strain and the two reductions by the
// same Lee-Fenves weight r that already splits the two hardening measures, so
// the two bookkeepings cannot disagree, and the two accumulators sum to
// lambda*m exactly. The dilation angle survives the split EXACTLY: sharing the
// multiplier scales m uniformly and psi is the RATIO of its volumetric to its
// deviatoric part. Only the AMOUNT of flow is shared, never its direction.
//
// r IS READ ON THE EFFECTIVE STRESS AND NOT ON THE NOMINAL ONE, which is the
// single most expensive thing to get wrong here. r splits a rate of the
// EFFECTIVE space, so the two reductions cannot be inside their own splitting
// rule; read on the nominal stress, a 600-fold asymmetry between omega_t and
// omega_c on a cracked tensile leg promotes a 0.26 MPa transverse compression
// over a 45 MPa tension, the model calls a purely tensile state 97% compressive,
// makes 68% of the flow permanent on a leg whose dial promises none, advances
// kappa_c in tension, and the return mapping limit-cycles.
//
// UNDER IMPL-EX ONE SCALAR IS EXTRAPOLATED - the plastic multiplier - and
// everything else follows from it. The spectral split PT/PC is FROZEN, the two
// reductions are NOT: a projector is not a monotone function of an irreversible
// measure and an extrapolated one means nothing, while omega is a function of
// kappa, which IS the measure this scheme extrapolates. Costing nothing:
// dlambda comes from the COMMITTED state, so kt and kc are constants of the
// step, sigma stays affine in eps and the algorithmic tangent is exactly W:C.
//
// THE MATERIAL MEASURES THE IMPL-EX ERROR AND NEVER DECIDES. The metric is the
// normalized stress gap |s_implex - s_implicit|/s_ref, the same quantity the
// whole family measures, with the largest Voigt component as the norm and the
// strongest stress the material can carry as the denominator. Rejecting a step
// belongs to CTestImplexWrapper: a material that aborts owns a policy it cannot
// see the analysis to choose.
//
// VOIGT CONVENTION, and it is the decision that shapes half the implementation.
// Every 6-vector INSIDE this material carries the TENSOR components of the shear
// (eps_12, not gamma_12 = 2 eps_12). The engineering convention appears in
// exactly two functions - setTrialStrain, which halves the last three
// components, and the two tangents, which halve the last three COLUMNS. The
// reason is that everything the model does - the eigen-decomposition, the
// Macaulay brackets on principal stresses, the deviator, the flow direction - is
// natural on tensor components and full of factor-of-two traps otherwise, and
// there are twenty such places against two boundaries. Two things fall out for
// free: the elastic law is a closed form, lam*tr(e)*I + 2*mu*e, with no matrix
// at all; and a fourth-order operator becomes a 6x6 whose shear COLUMNS carry
// the factor two of the contraction over the off-diagonal pairs, which is
// exactly what ASDSpectralSplit::computePjj already builds.
//
// WHAT THIS MODEL DOES NOT DO, stated so nobody looks for it:
//
// * no nonlinear unloading, by request. Both unloadings are straight lines; what
//   the dial chooses is WHERE each one lands. A side at df = 1 has no hysteresis
//   at all, which is what "pure damage" means;
// * the unilateral defect survives, but BOUNDED: a compressive eps_p makes sbar
//   positive at eps = 0, so a crushed specimen unloaded past its plastic strain
//   goes into tension - and that spurious tension is a positive eigenvalue, so
//   it is capped by the tensile mechanism at the tensile backbone instead of
//   running to E*|eps_p|;
// * THE TANGENT IS NOT SYMMETRIC - 24% of E in uniaxial tension - and not
//   positive definite everywhere: on the approach to the tensile vertex the
//   symmetric part reaches -1.4e5. That is the formulation and not the
//   integration (it predates the line search by the same numbers), and a
//   softening tangent losing definiteness is not by itself pathological. Under
//   IMPL-EX the tangent is exact and none of this arises;
// * THE TENSILE APEX is where the local gradient stops being usable: a dilatant
//   non-associated flow pushes the stress onto the hydrostatic tensile axis,
//   where the deviator vanishes, m degenerates to tan(psi)/3*I and the surface
//   closes at qc(1-alpha)/(3 alpha + beta) - 0.31 MPa with a calibrated
//   concrete, i.e. effectively "no hydrostatic tension". There the consistency
//   equation is still a scalar equation on a bounded interval, so it is
//   BRACKETED instead, which needs a sign change and no derivative at all. At
//   psi = 0 even that is unavailable and is reported as a failed step, honestly:
//   a purely deviatoric flow has m identically zero on the hydrostatic axis, so
//   no multiplier changes the stress and the residual is not a function of
//   lambda.
//
// NO NUMERICAL TANGENT FLAG, unlike ASDConcrete3D, and it is not an omission.
// Uniaxial compression sits on TWO non-smooth features of the surface at once -
// a doubly degenerate s_max and the kink of the Macaulay pair - so a one-sided
// difference resolves the degeneracy in whichever direction it happens to
// perturb and reports a derivative off by a factor of two in a dominant term.
// Measured on the reference implementation: with the numerical tangent in place
// the mixed-control Newton oscillated forever and 181 of 300 steps of the
// uniaxial compression path failed to converge while the return mapping ITSELF
// was converging in two iterations at every one of them. The analytic form
// carries the averaged subgradient through, so the Jacobian sees the same
// symmetric state the material does.
//
// The reference implementation, and the oracle this was ported against, is
// openseesideas/asd_cdp_damage_3d.py (with its base asd_cdp_3d.py) in the
// OpenSeesIDEAS bench.

#ifndef ASDPlasticDamageConcrete3DMaterial_h
#define ASDPlasticDamageConcrete3DMaterial_h

#include <NDMaterial.h>
#include <ASDHardeningLaw.h>
#include <IMPLEXManager.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>
#include <cmath>

/**
q(kappa) and dq/dkappa: a 1D backbone read as a CDP hardening curve.

Built from an ASDHardeningLaw by the map kappa = x - q/E, which is monotone
exactly under the admissibility conditions ASDHardeningLaw::adjust() already
imposes - non-decreasing plastic strain, tangent never above E. So the
conditions this family has always imposed on its backbones ARE the conditions
for them to be legal CDP hardening curves. The origin is dropped: the curve
starts at kappa = 0 carrying the ELASTIC LIMIT, which is where a yield surface
starts.

THE DAMAGE CAP IS 0 IN THIS MODEL, which is what makes the curve the NOMINAL
backbone: abscissa x - y/E, ordinate y, and the input damage arrays ignored by
construction. The cap is a parameter and not a constant because the classical
CDP core needs 1 (or 0.99: uncapped, the Model Code tensile array ends at
d = 0.999993 and asks for q = 434.85 MPa on a 3 MPa material, with two adjacent
points a hardening modulus of 4e20 apart).

Points that do not advance kappa are DROPPED and counted. They are not
interpolable - q jumping at constant kappa is an infinite hardening modulus -
and silently smoothing them would hide the one thing a user needs to be told
about their input.
*/
class ASDCDPHardeningCurve
{
public:
	ASDCDPHardeningCurve() = default;

	// Builds the curve. Returns false, and leaves a reason in 'why', when the
	// law cannot be used: no usable point, or a non-positive effective strength
	// at any of them. THE SECOND ONE IS WORTH CHECKING HERE rather than
	// downstream: beta = qc/qt*(1-a) - (1+a) has qt in a denominator and dF/dkt
	// has qt*qt, so a zero would surface as a division by zero from inside the
	// return mapping.
	bool build(const ASDHardeningLaw& law, double E, double damage_cap,
		const char** why);

	// (q, dq/dkappa), linear between points.
	//
	// Beyond the last point the value is held CONSTANT - perfectly plastic in
	// effective stress - rather than extrapolated along the last tangent, which
	// is what the 1D law does for a positive one. Holding is the only safe
	// choice in both directions here, because an extrapolated negative tangent
	// would drive the effective strength through zero and take the yield
	// surface with it.
	//
	// AT the first point the slope is the FORWARD one and not zero: a plastic
	// step starting from kappa = 0 has to see the hardening it is about to walk
	// into. Returning zero there made the first corrector of every return
	// mapping miss the whole hardening term of its denominator - measured, it
	// overshot the compressive yield by 0.29 MPa and then needed a NEGATIVE
	// corrector to come back.
	void evaluate(double kappa, double& q, double& dq) const;

	// The NOMINAL strength at kappa, against evaluate's effective one. Equal to
	// it when the backbone carries no damage, which is this model's case.
	double nominal(double kappa) const;

	// y/q = 1 - d at kappa: the reduction the ENVELOPE carries, as opposed to
	// the one this model derives from the dial. 1 for an undamaged backbone.
	double omegaOfCurve(double kappa) const;

	inline bool isValid() const { return m_kappa.size() > 0; }
	inline double initialYield() const { return m_q.size() ? m_q[0] : 0.0; }
	inline double maxEffectiveStress() const { return m_qmax; }
	inline double lastKappa() const { return m_kappa.size() ? m_kappa.back() : 0.0; }
	inline std::size_t size() const { return m_kappa.size(); }
	inline int dropped() const { return m_dropped; }
	inline const std::vector<double>& kappa() const { return m_kappa; }
	inline const std::vector<double>& q() const { return m_q; }
	inline const std::vector<double>& y() const { return m_y; }

private:
	std::vector<double> m_kappa;
	std::vector<double> m_q;
	std::vector<double> m_y;
	double m_qmax = 0.0;
	int m_dropped = 0;
};

class ASDPlasticDamageConcrete3DMaterial : public NDMaterial, public IMPLEXObject
{
public:
	using HardeningLawPoint = ASDHardeningLawPoint;
	using HardeningLawType = ASDHardeningLawType;
	using HardeningLawPointComponent = ASDHardeningLawPointComponent;
	using HardeningLaw = ASDHardeningLaw;

	/**
	Everything the implicit pass of integrate() can write. Saved and restored
	around the non-destructive measurement of the IMPL-EX error, so that a step
	that is measured and then rejected is a step that never happened.

	IT HAS TO BE COMPLETE, and revertToLastCommit() is not a substitute: it does
	not restore the FROZEN spectral split, which is a committed quantity that an
	implicit pass overwrites. Two fields beyond what the reference
	implementation's peek carries - the integration diagnostics and C - because
	the recorders run between convergence and commit and must see what the step
	delivered, not what the measurement computed.
	*/
	struct TrialState {
		Vector stress = Vector(6);
		Vector sbar = Vector(6);
		Vector sbar_pos = Vector(6);
		Vector sbar_neg = Vector(6);
		Vector ep_cr = Vector(6);
		Vector ep_pl = Vector(6);
		Vector d_split = Vector(3);
		Matrix V_split = Matrix(3, 3);
		Matrix W = Matrix(6, 6);
		Matrix C = Matrix(6, 6);
		double kt = 0.0;
		double kc = 0.0;
		double dlambda = 0.0;
		double implex_error = 0.0;
		int n_iter = 0;
		int n_residual = 0;
		int bisected = 0;
		bool plastic = false;
		bool failed = false;
		bool stagnated = false;
		bool extrapolated = false;
	};

public:
	// life-cycle
	ASDPlasticDamageConcrete3DMaterial(
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
		double _stagnation_tol);
	ASDPlasticDamageConcrete3DMaterial();
	~ASDPlasticDamageConcrete3DMaterial();

	// info
	const char* getClassType(void) const { return "ASDPlasticDamageConcrete3DMaterial"; };

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
	virtual int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);

	// IMPL-EX error control
	double computeImplexErrorMetric(void);
	double implexTimeRatio(void) const;

private:
	// --- setup ------------------------------------------------------------ //

	// Regularize once, then convert the two backbones into CDP hardening
	// curves. LAZY, and it has to be: regularization needs the parent element's
	// characteristic length, which is only known at the first setTrialStrain.
	// So a hardening table that cannot be converted can only be reported from
	// there - which is why the parser does a dry run of the same conversion on
	// the un-regularized laws, where the user can still act on it.
	bool prepare(void);

	// --- the surface and the potential ------------------------------------ //

	// F(s, qt, qc), in stress units, negative inside. With TENSION POSITIVE and
	// I1 = tr(s); Abaqus writes it with the equivalent pressure p = -I1/3, so
	// its -3*alpha*p is this alpha*I1:
	//
	//   F = 1/(1-a) * ( q + a*I1 + beta*<s_max> - gamma*<-s_max> ) - qc
	//
	// It yields at sc = qc in uniaxial compression, at st = qt in uniaxial
	// tension - that is what beta is FOR, and why the surface needs both
	// hardening variables at once - and at fb0/fc0 times the uniaxial strength
	// in equibiaxial compression, which is what alpha was solved from.
	double yieldFunction(const Vector& s, double qt, double qc) const;

	// qc/qt*(1-a) - (1+a), and the ONE place qt is divided by.
	//
	// THE FLOOR IS ON REPRESENTABILITY, NOT ON STRENGTH. A floor inside the
	// tabulated range - qt at ft/10, which is what an earlier cure did -
	// truncates the softening tail and leaves a residual tensile strength the
	// input never asked for. m_q_tiny is 1e-24 times the largest tabulated
	// strength: twenty-four orders of magnitude below anything the curve will
	// accept, so it cannot alter a value the model computes. What it does is
	// keep beta FINITE, and that matters because beta multiplies the
	// convergence tolerance through residualScale(): an inf there makes
	// |F| <= tol*scale true for every state, and the return mapping then
	// reports a converged step it never took. A silent false convergence is the
	// worst failure mode available and it costs one max to make it unreachable.
	double surfaceBeta(double qt, double qc) const;

	// (dF/dkt, dF/dkc). BOTH effective strengths move, and so does beta, which
	// depends on both - that coupling is what "two hardening variables on one
	// surface" means.
	void dFdKappa(double smax, double qt, double qc, double dqt, double dqc,
		double& dF_dkt, double& dF_dkc) const;

	// dF/ds, with the two non-smooth places handled explicitly. It hands back
	// s_max as well, because it has just paid for the decomposition that
	// produces it and the caller needs exactly that value for dFdKappa: the
	// alternative is a second eigen-decomposition of the same stress per
	// iteration, which is what the reference implementation does.
	void dFds(const Vector& s, double qt, double qc, Vector& out, double& smax) const;

	// (s_max, d s_max / ds), AVERAGED over a repeated maximum.
	//
	// This is the sharpest edge in the whole formulation and it sits exactly on
	// the load case the model is calibrated against. d s_max / ds is n (x) n
	// only when the largest principal stress is SIMPLE; when it is repeated the
	// derivative does not exist and its subdifferential is the convex hull of
	// the projectors of the whole eigenspace. Uniaxial compression has
	// s = (-sc, 0, 0): the maximum is 0 with multiplicity two, so an
	// eigen-solver returns an ARBITRARY direction of the 22/33 plane and the
	// gradient comes out asymmetric in two directions the load case cannot
	// distinguish. Measured before the averaging: the return mapping failed at
	// 228 of 300 steps past the elastic limit in uniaxial compression and the
	// stress ran 211 MPa off the backbone.
	void smaxProjector(const Vector& s, double& smax, Vector& P) const;

	// dq/ds = 3S/(2q), computed so that no stress SCALE can lose it. The
	// quantity is scale-INVARIANT - sqrt(3/2) times the unit deviator - so
	// normalizing before the ratio changes no value and removes the only way
	// the expression can fail: the direct form squares S inside q, so a
	// deviator of 1e-170 underflows the sum to zero, the guard reads q == 0 and
	// reports NO deviatoric gradient for a state that has a perfectly good one.
	// At S = 0 exactly the Mises cone genuinely has no gradient - the
	// hydrostatic axis is the only place it happens - and zero is the symmetric
	// element of its subdifferential.
	void dqds(const Vector& s, Vector& out) const;

	// m = dG/ds of the hyperbolic Drucker-Prager potential,
	//
	//   G = sqrt( (ecc*ft0*tan(psi))^2 + q^2 ) + I1/3 * tan(psi)
	//
	// with tr(m) = tan(psi) EXACTLY at every stress state: the dilation angle is
	// the volumetric plastic strain per unit deviatoric flow, with nothing in
	// between to calibrate. The hyperbola is not decoration: it is what gives
	// the surface a well-defined gradient at the tensile apex, where an
	// associated Drucker-Prager has a cone point.
	void flowDirection(const Vector& s, Vector& out) const;

	// Lee and Fenves' r: HOW TENSILE the state is, in [0, 1], read ON THE
	// EFFECTIVE STRESS - see the note in the file header on what reading it on
	// the nominal one costs.
	double splitWeight(void) const;

	// (h_t, h_c): the advance of each measure per unit dlambda, by Lee and
	// Fenves' split. The Macaulay brackets are what keeps both measures
	// non-decreasing, which is a thermodynamic requirement and not a guard.
	void hardeningRates(const Vector& m, double r, double& h_t, double& h_c) const;

	// --- the reductions --------------------------------------------------- //

	// q / (df*E*kappa + q): the reduction a side has earned. DERIVED, for every
	// df and not only at the two ends - see the file header.
	double omega(double kappa, double df, const ASDCDPHardeningCurve& hard) const;
	// d omega / d kappa = df*E*(dq*kappa - q) / z^2. Negative wherever the
	// secant modulus of the backbone decreases, which is what makes the damage
	// terms of the denominator help rather than fight.
	double domega(double kappa, double df, const ASDCDPHardeningCurve& hard) const;
	// how much of dlambda*m becomes PERMANENT strain: (1-df_t) of the tensile
	// share and (1-df_c) of the compressive one, weighted by r. The rest is
	// carried by the two reductions instead, and the two sum to dlambda*m
	// exactly.
	inline double plasticShare(double r) const {
		return (1.0 - damage_t) * r + (1.0 - damage_c) * (1.0 - r);
	}

	// W from a FROZEN spectral record and the two measures. Two callers - the
	// explicit IMPL-EX pass and revertToLastCommit - and they must not be able
	// to disagree about the operator they rebuild.
	void rebuildOperator(const Vector& d_rec, const Matrix& V_rec,
		double kt_, double kc_);

	// The damaged elastic law, and the place the split is built: writes sbar,
	// the split record, W, sbar_pos, sbar_neg, and returns the NOMINAL stress
	// in 'out'. ep_pl is the plastic strain; ep_cr, the cracking accumulator,
	// is NOT part of the elastic law at all and is kept for output.
	int effective(const Vector& eps, const Vector& ep_pl, double kt, double kc,
		Vector& out);

	// --- the integration -------------------------------------------------- //

	// What |F| has to be measured against, and it is not qc. F is not a
	// distance to the surface: it is a distance times the local gradient, and on
	// the tensile meridian that gradient is beta/(1-alpha) ~ qc/qt. With a
	// softening tensile curve reaching ft/1000 the amplification is FOUR
	// THOUSAND, so an absolute test |F| <= tol*qc asks for the stress to 3e-13
	// MPa while the effective stress carrying it is 41 MPa - below double
	// precision, unreachable, and the return mapping stalls at max_iter with a
	// converged multiplier it is not allowed to accept. Dividing by the
	// amplification puts the tolerance back in stress units.
	double residualScale(void) const;

	// |F| converted to a DISTANCE TO THE SURFACE, in stress units.
	//
	// F is a stress distance TIMES the local gradient, and on the tensile
	// meridian that gradient is 1 + |beta|/(1-alpha) ~ qc/qt. Dividing it out
	// gives the only number a VERDICT can honestly be made on: how far outside
	// its own surface the state about to be handed back actually is.
	//
	// WHERE THIS IS USED AND WHERE IT IS NOT, because the distinction is the
	// whole design. It decides the verdict - failed against stagnated - and it
	// does NOT decide acceptance. The iteration still stops exactly where it
	// stopped before, so every step that converges converges to the same bits.
	// Measured on the bench over seven dial pairs and three combinations: the
	// delivered stress and both hardening variables move by EXACTLY zero, while
	// the failure counts fall (58 -> 9 at damage 0.3/0, 4 -> 0 at 0.5/0.3,
	// 1 -> 0 at 1/0.3). Only the flag moves, and the flag is what a host acts
	// on: a -1 out of setTrialStrain makes the analysis cut its step.
	//
	// WHY IT WAS NEEDED, from the trace of the one step at (1, 0.3) that used to
	// reach max_iter. Purely tensile on the softening tail, qt = 0.0306,
	// qc = 27.17, so the amplification is 886.9. The iteration neither diverges
	// nor wanders - by the twentieth pass the frozen rates have settled and F
	// sits on a fixed value - and the best iterate it saw was F = 0.0201, i.e.
	// 2.3e-5 MPa from the surface. The acceptance test asks |F| <= 3e-9, which
	// at that amplification is 3.4e-12 MPa: below what double precision can
	// express next to a 40 MPa effective stress. The step had a usable answer
	// and called itself failed.
	double residualStress(double f, double qt, double qc) const;

	// A = -d sigma / d lambda: one plastic term and TWO damage terms. The
	// damage terms are the ones a classical CDP has no place for; they move the
	// stress with NO strain, which is precisely how a side whose df is 1 is
	// brought back to the surface without leaving anything permanent, and at
	// df = 0 the corresponding domega is identically zero and the term
	// disappears. Dropping one is a first-order error in the denominator: on the
	// model this was ported from, leaving such a term out took the return
	// mapping from 4 iterations to 40, with failures.
	void corrector(double r, const Vector& m, double h_t, double h_c,
		double kt, double kc, Vector& out) const;

	// The largest multiplier this step is allowed to try. Not a convergence
	// device: a BRACKET, so the line search and the bisection have a finite
	// interval and a runaway iterate is caught where it starts instead of after
	// it has written kt = 1.88e9 into the state.
	double lambdaCap(double h_t, double h_c, const Vector& m,
		const Vector& eps) const;

	// Is this state the TENSILE VERTEX, rather than merely a hard one? den <= 0
	// does not answer it, and the difference decides whether the bracketing
	// fallback is allowed to move the state a long way. Measured over the
	// alternating protocol, all 49 den <= 0 events were shear-dominated - q
	// between 32 and 67 MPa against a mean stress of 1e-3 - and not one was near
	// a vertex; handing those to the bracket let it jump two decades in lambda,
	// and the distant root of the FROZEN equation is not a root of the real one.
	// The vertex is where the deviator has nothing left to say: mean stress
	// tensile, and the Mises equivalent no larger than it.
	bool atApex(const Vector& s) const;

	// The state a multiplier produces, and its residual. ONE place, so that the
	// line search, the bisection and the iterate finally adopted cannot disagree
	// about what lambda means. It leaves W / sbar_pos / sbar_neg holding the
	// operator of the state it just built, which is what corrector() and
	// getTangent() read - hence the rule, enforced by the loop, that the
	// ADOPTED trial must be the LAST ONE EVALUATED.
	struct Trial {
		Vector s = Vector(6);
		Vector ep_cr = Vector(6);
		Vector ep_pl = Vector(6);
		double kt = 0.0;
		double kc = 0.0;
		double f = 0.0;
		double ftol = 0.0;
		bool ok = true;
	};
	void trialAt(double lam, const Vector& eps, const Vector& m, double pf,
		double h_t, double h_c, Trial& out);

	// Bracket the consistency equation on (lam0, cap] and bisect it: the
	// fallback for the tensile apex and for nothing else. It gives up exactly
	// one thing, the quadratic rate, and in exchange needs no derivative at all,
	// which is the only honest position to take at a vertex of the surface. The
	// scan is GEOMETRIC and upward, from cap*2^-20 to cap, because the
	// multiplier a step needs spans the whole range: the elastic-limit steps of
	// a fine ramp want 1e-9 and an apex event wants the largest admissible
	// value. Returns false when the residual keeps its sign over the whole
	// interval, which is the real answer for a hydrostatic state at psi = 0.
	bool bracket(double f0, double lam0, double cap, const Vector& eps,
		const Vector& m, double pf, double h_t, double h_c,
		double& lam_out, Trial& out);

	// The return mapping. Cutting plane with the damage terms in the
	// denominator, a damped corrector and the bracket as a fallback; the update
	// is BACKWARD EULER (total multiplier, direction at the current iterate) and
	// not an accumulation of the iteration path - see the .cpp.
	int integrate(const Vector& eps);

	// --- the tangent ------------------------------------------------------ //

	// Build C from the state the last pass left. Called at the end of
	// setTrialStrain, so that getTangent() is a plain accessor as it is
	// everywhere else in the family, and so that the operator handed to the
	// element belongs to the state the element was given - which under
	// IMPL-EX with the in-material control means the DELIVERED one, after the
	// peek has undone its implicit pass.
	void computeTangent(void);

	// --- IMPL-EX ---------------------------------------------------------- //

	// dt_n / dt_n_commit * alpha, CLAMPED AT ZERO. The clamp is not cosmetic: a
	// negative ratio extrapolates an irreversible process BACKWARDS, which is
	// not a state this material can be in - the plastic strain stops being
	// monotone along the flow direction and the hardening measures decrease,
	// undoing damage. The precondition it guards is that dt measures the SIZE
	// of the imposed increment, which is stronger than "positive": halve the
	// increment without halving dt and the extrapolation over-predicts by two,
	// silently, with the error metric blaming the discretization for it.
	double timeFactor(void) const;

	// The explicit pass: one scalar extrapolated, the split FROZEN, the two
	// reductions NOT - see the file header.
	void extrapolate(const Vector& eps);

	// |s_implex - s_implicit|/s_ref, the family's metric.
	double implexStressGap(const Vector& delivered, const Vector& stress_implicit) const;
	// the largest stress this material can carry: the metric's denominator, so
	// that one tolerance means the same thing at every load level and in every
	// model of the family
	double stressReference(void) const;

	void saveTrialState(TrialState& x) const;
	void restoreTrialState(const TrialState& x);

	// --- responses -------------------------------------------------------- //
	Vector getHardeningLawVector(HardeningLawType ltype, HardeningLawPointComponent c) const;
	Vector getCDPCurveVector(HardeningLawType ltype, int component) const;
	const Vector& getKappa() const;
	const Vector& getDamage() const;
	const Vector& getOmega() const;
	const Vector& getStrength() const;
	const Vector& getSplitWeight() const;
	const Vector& getPlasticStrainVector() const;
	const Vector& getCrackingStrainVector() const;
	const Vector& getEffectiveStress() const;
	const Vector& getImplexError() const;
	const Vector& getImplexStress() const;
	const Vector& getTimeIncrements() const;
	const Vector& getIntegrationInfo() const;

private:
	// --- material parameters ---------------------------------------------- //
	// Young's modulus
	double E = 0.0;
	// Poisson's ratio
	double nu = 0.2;
	// mass density
	double rho = 0.0;
	// the dilation angle, in RADIANS in the meridional plane at high confining
	// pressure. The user parameter this whole change of formulation is for
	double psi = 0.0;
	// eccentricity of the hyperbolic potential. 0 is the sharp Drucker-Prager
	// cone, 0.1 is Abaqus' default
	double ecc = 0.1;
	// ratio of the equibiaxial to the uniaxial compressive strength. alpha is
	// solved from it
	double fb0_fc0 = 1.16;
	// ratio of the second stress invariant on the tensile meridian to that on
	// the compressive one. 2/3 is Abaqus' default, 1 makes the deviatoric
	// section a circle
	double Kc = 2.0 / 3.0;
	// Lubliner's alpha, from fb0_fc0
	double alpha = 0.0;
	// the third-invariant shape factor 3(1-Kc)/(2Kc-1), active only where
	// s_max < 0
	double gam = 0.0;
	// the two dials, one per side: 0 = pure plasticity, 1 = pure damage
	double damage_t = 1.0;
	double damage_c = 0.0;

	// --- IMPL-EX ---------------------------------------------------------- //
	// True = use the IMPL-EX algorithm
	bool implex = false;
	// True = measure the error at every setTrialStrain and publish it. The
	// LEGACY path: it pays one implicit solve per Newton iteration, at strains
	// that are not equilibrated yet. The measurement the error control runs on
	// is the one CTestImplexWrapper asks for, once per step, through
	// computeImplexErrorMetric()
	bool implex_control = false;
	// True = let this material fail the step on its own when the error is over
	// tolerance. OFF by default and it should stay off: a material that aborts
	// owns a policy it cannot see the analysis to choose
	bool implex_abort_on_error = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	double implex_alpha = 1.0;

	// --- integration controls --------------------------------------------- //
	// relative tolerance of the return mapping, against residualScale()
	double tol = 1.0e-10;
	int max_iter = 100;
	// how many times the Newton corrector may be halved before the step is
	// handed to the bracketing fallback
	int max_backtrack = 12;
	// the residual, relative to residualScale(), at which an iteration that can
	// no longer reduce |F| is ACCEPTED instead of reported as a failure. F is
	// only continuous across a change of sign of an eigenvalue, so tol is not
	// always reachable and the honest thing is to say where the iteration
	// stopped. 0 restores the all-or-nothing behaviour
	double stagnation_tol = 1.0e-6;
	// How far above the residual the STEP STARTED FROM a line-search trial may
	// sit and still be accepted. 1 reads "a corrector may not leave the state
	// worse than it found it", which is much weaker than a monotone descent test
	// and is the whole safeguard. Not a material property, hence not an input
	static const double ResidualGrowth;

	// --- the backbones ---------------------------------------------------- //
	bool auto_regularize = false;
	bool regularization_done = false;
	// the CDP curves exist. A SEPARATE flag from regularization_done, and it has
	// to be: on the receive path the laws arrive already regularized, so the
	// regularization must be skipped while the conversion must still happen
	bool curves_ready = false;
	bool setup_ok = false;
	bool setup_reported = false;
	double lch = 1.0;
	double lch_ref = 1.0;
	HardeningLaw ht;
	HardeningLaw hc;
	// the two backbones as CDP hardening curves. Built lazily by prepare(),
	// because regularization needs the element's characteristic length
	ASDCDPHardeningCurve hard_t;
	ASDCDPHardeningCurve hard_c;
	// REPRESENTABILITY floor on qt - see surfaceBeta(). The starting value looks
	// over-large for a floor and is the smallest safe one: dF/dkt divides by
	// qt*qt, so a floor below 1.5e-154 squares to zero and reintroduces exactly
	// the division it is there to prevent
	double q_tiny = 1.0e-150;

	// --- state ------------------------------------------------------------ //
	// TENSOR components of the shear throughout - see the file header
	Vector strain = Vector(6);
	Vector strain_commit = Vector(6);
	// what getStrain() hands back, in the ENGINEERING convention. A member so
	// that an element may hold the reference across integration points - see
	// getStrain(). Not state: it is refreshed on every call
	Vector strain_eng = Vector(6);
	// the NOMINAL stress: what equilibrium sees and what the surface is written
	// in
	Vector stress = Vector(6);
	Vector stress_commit = Vector(6);
	// WHAT THE STEP DELIVERED TO THE ELEMENT, kept because it is otherwise
	// unobservable: under IMPL-EX commitState re-solves implicitly and installs
	// that answer over `stress`, so from outside - a recorder, a test - the
	// extrapolated stress is gone by the time anyone can ask. It is what the
	// error metric measures the distance from, and publishing it is what lets
	// that distance be checked component by component instead of in norm. On the
	// implicit path there is no extrapolation and it equals `stress`
	Vector stress_implex = Vector(6);
	// the effective stress. Kept alongside the nominal one because splitWeight()
	// answers from it and the explicit pass never calls effective()
	Vector sbar = Vector(6);
	Vector sbar_pos = Vector(6);
	Vector sbar_neg = Vector(6);
	// the cracking accumulator: the share of the flow the REDUCTIONS carry.
	// With the plastic strain it sums to lambda*m exactly, which is the
	// bookkeeping invariant. Not part of the elastic law
	Vector ep_cr = Vector(6);
	Vector ep_cr_commit = Vector(6);
	// the plastic strain: the share of the flow that is permanent
	Vector ep_pl = Vector(6);
	Vector ep_pl_commit = Vector(6);
	// the two hardening measures
	double kt = 0.0;
	double kc = 0.0;
	double kt_commit = 0.0;
	double kc_commit = 0.0;
	// the damaged operator W = omega_t*PT + omega_c*PC, of the state that was
	// last built. ONE member and not two: the reference implementation keeps a
	// per-iteration W alongside the adopted one, and the loop's own invariant -
	// the adopted trial is the last one evaluated - makes them equal everywhere
	// either is read
	Matrix W = Matrix(6, 6);
	// the tangent, as getTangent() must hand it back
	Matrix C = Matrix(6, 6);
	// THE SPECTRAL SPLIT IS NOT STORED AS TWO 6x6 MATRICES, it is stored as the
	// decomposition it comes from: 12 doubles instead of 72, per Gauss point,
	// and the rebuild is the same function of the same data so it cannot drift.
	// d_split holds the eigenvalues (the sign band needs them, not just the
	// signs) and V_split the directions in columns
	Vector d_split = Vector(3);
	Matrix V_split = Matrix(3, 3);
	// the FROZEN decomposition the explicit pass rebuilds W on
	Vector d_split_commit = Vector(3);
	Matrix V_split_commit = Matrix(3, 3);
	// what the explicit pass extrapolates, and what it freezes. The direction is
	// the one the return path actually TOOK over the committed step, not the
	// flow direction at its end: they differ by the curvature of the return
	// path, and the first one is what makes the explicit update reproduce the
	// implicit one exactly when the extrapolated multiplier happens to be right
	double dlambda = 0.0;
	double dlambda_commit = 0.0;
	Vector m_commit = Vector(6);    // the plastic-strain rate
	Vector mcr_commit = Vector(6);  // the damage-carried strain rate
	double kt_rate_commit = 0.0;
	double kc_rate_commit = 0.0;
	double r_commit = 0.0;
	// implex time bookkeeping
	double dtime_n = 0.0;
	double dtime_n_commit = 0.0;
	double dtime_0 = 0.0;
	bool dtime_is_user_defined = false;
	bool commit_done = false;
	double implex_error = 0.0;
	// diagnostics of the last integration. n_residual is the number of yield
	// function evaluations, which is what the step actually costs; stagnated
	// means the residual stopped improving above tol but below stagnation_tol;
	// bisected counts the iterations that had to fall back on bracketing
	int n_iter = 0;
	int n_residual = 0;
	int bisected = 0;
	bool plastic = false;
	bool failed = false;
	bool stagnated = false;
	// True when the last setTrialStrain took the extrapolated path
	bool extrapolated = false;
};

#endif // ASDPlasticDamageConcrete3DMaterial_h
