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
// A damage-free hysteretic material with a shapeable unload/reload law.
//
// PEAK-ORIENTED, AND WITHOUT A DAMAGE VARIABLE. In the plasticity-based
// formulation of ASDConcrete1D the unilateral effect - a stiffness that comes
// back when the crack closes - could only be had by carrying two damage
// variables, one per side: the (1-d)E secant was the only device able to bring
// the response back toward the opposite backbone. The peak-oriented rule makes
// that device unnecessary, because the reloading branch aims at the EXACT point
// of the opposite backbone whatever permanent strain was left behind, so the
// unilateral effect is produced by the branch topology instead. With d = 0
// everywhere the effective and the nominal backbone stress coincide, the
// permanent strain of a backbone point is x - y/E, and the Td/Cd damage arrays
// are refused rather than ignored: a non-zero damage would silently reshape the
// backbone inside ASDHardeningLaw::adjust().
//
// THE SHAPE OF THE UNLOAD/RELOAD PATH. The excursion that leaves a point S of
// one side's backbone and ends on the opposite target T is the polyline
//
//     S  ->  P1 = (c1, 0)  ->  V_r  ->  T
//
// with c1 the side's own plastic strain (the descent is the elastic line) and
// V_r the shaping vertex of the reloading half. a4 is the height of the vertex
// in the reloading branch's stress span, a3 its horizontal position at that
// height, sliding between three references ordered from the stiffest first leg
// to the softest: the prolongation of the last unloading leg (a3 = -1), the
// peak-oriented chord (a3 = 0, the previous law), the elastic line through the
// target (a3 = +1). Two HARD CONSTRAINTS are enforced whatever the parameters:
// the first reloading leg is never steeper than the last unloading one, and the
// vertex never passes the elastic line through the target. So every leg has a
// slope in (0, E], the operator is always positive definite and never stiffer
// than the elastic one, and the path still ends exactly on T.
//
// THREE CORNERS, THREE CUTS. The corner at P1 - where the unloading of one side
// hands over to the reloading of the other - is cut by the chamfer, and that is
// what gives permanent strain back. The reloading vertex and the junction where
// the reload runs back into the envelope are cut too (smooth, tip_in/tip_out).
// Every offset is a FRACTION OF ITS OWN LEG, so no length in the sigma-eps plane
// - and therefore no choice of units - ever has to be defined, and a fraction
// <= 1 can never overrun the leg it is measured on.
//
// NOTHING PERSISTS. The branches are built inside compute() and dropped when it
// returns; the state carried between steps is scalars only, plus the two
// backbones and the two shapes. So the whole geometry machinery adds not one
// byte per Gauss point, which is what matters where a copy of this material
// lives on every fibre of every section of every integration point. The branch
// is a FIXED-CAPACITY buffer for the same reason: the reference implementation
// allocates one list per branch per Gauss point per iteration, and it does not
// need to - the number of points is bounded, see ASDHysteretic1DBranch.
//
// IMPL-EX. The equivalent strain measures xt/xc are extrapolated from the last
// two committed states; every on/off state is FROZEN at its last committed
// (implicitly computed) value: the active side, the branch types, their birth
// points, the measure a reload aims at, and WHICH LEG of the branch is active -
// the last one without a new state variable, because the leg is the one holding
// the committed strain. Within a step the geometry and the leg are therefore
// constant, sigma is linear in the strain and the algorithmic tangent is the
// leg's slope: step-wise linearity is preserved exactly, whatever the
// tessellation of a rounded corner.
//
// THE ERROR METRIC IS THE NORMALIZED STRESS GAP, |s_implex - s_implicit|/s_ref,
// the same quantity ASDConcrete1D and ASDConcrete3D now measure, so that one
// tolerance means the same thing across the family. The material MEASURES and
// never decides: rejecting a step belongs to CTestImplexWrapper.
//
// USER LIMIT STATES. '-limitStates $x1 $x2 ...' takes an ascending list of
// abscissae of the strain measure - IO/LS/CP, or however many the user wants -
// and the 'limitStateRatio' response reports the HIGHEST one ever reached, as a
// continuous position in that scale: 0 = none, 1 = exactly the first, 2.5 =
// halfway between the second and the third, n = the last. ONE scalar, and in the
// LS scale rather than in strain, because the i-th limit state sits at a
// different abscissa on each side and in every other section of the model while
// its index does not - so the index is what can be compared and maximized across
// them. It costs no state variable: xt/xc are monotone and already ARE the peak
// reached on their side, which also makes the answer follow a rollback for free.
// See getLimitStateRatio(). '-limitStatesC' overrides the list on the
// compressive side, for an asymmetric backbone.
//
// The reference implementation, and the oracle this was ported against, is
// openseesideas/asd_hysteretic_1d.py in the OpenSeesIDEAS bench.

#ifndef ASDHysteretic1DMaterial_h
#define ASDHysteretic1DMaterial_h

#include <UniaxialMaterial.h>
#include <ASDHardeningLaw.h>
#include <IMPLEXManager.h>
#include <Vector.h>
#include <cmath>
#include <vector>

/**
Which sub-envelope branch a side is currently on.
*/
enum class ASDHysteretic1DBranchType {
	UnloadEnv = 0,   // unloading off the envelope, down to (+-c1, 0)
	UnloadPt = 1,    // unloading born at a reversal point, down to (+-c1, 0)
	ReloadCross = 2, // born at a zero crossing, aims at the side's target
	ReloadRev = 3    // born at a reversal, aims at the side's target
};

inline bool ASDHysteretic1DIsReload(ASDHysteretic1DBranchType t) {
	return t == ASDHysteretic1DBranchType::ReloadCross ||
		t == ASDHysteretic1DBranchType::ReloadRev;
}

/**
The shape of the unload/reload path of ONE side.

The nine fields are the GEOMETRY of the path. The way it is meant to be used is
one of the three named constructors, which derive them from parameters with a
mechanical meaning; all-zeros is the plain bilinear peak-oriented law, exactly
and branch by branch.
*/
class ASDHysteretic1DShape
{
public:
	// abscissa of the reloading vertex, in [-1, +1] (see the file header)
	double a3 = 0.0;
	// height of the reloading vertex in the branch's stress span, in [0, 1]
	double a4 = 0.0;
	// cut of the corner at zero stress, as a fraction of EACH of the two legs
	double chamfer = 0.0;
	// asymmetry of that cut, in [-1, +1]. Positive cuts the incoming
	// (unloading) leg further, so the path starts bending well above zero
	// stress instead of only at the corner - what a real concrete unloading
	// does
	double chamfer_bias = 0.0;
	// the same cut given DIRECTLY by its two ends instead, the way an
	// experimentalist reads them off a cyclic test: it starts at this fraction
	// of the stress the unloading started from ...
	double chamfer_rs = 0.0;
	// ... and ends this fraction of |DPE| past the corner, DPE being the
	// plastic strain an elastic unloading leaves. The only ABSOLUTE offset in
	// this class, so it is the only one that can ask for more leg than exists
	// and has to be capped
	double chamfer_dpe = 0.0;
	// round the cuts instead of chamfering them straight, in [0, 1]. 0 puts
	// the Bezier control point on the chord and the curve IS the straight cut;
	// 1 puts it on the corner and the curve is tangent to both legs
	double smooth = 0.0;
	// cut of the junction where the reload runs back into the envelope: how
	// far BEFORE the target it starts, as a fraction of the last reloading leg
	double tip_in = 0.0;
	// ... and how far PAST the target it ends, as a fraction of the same leg's
	// strain extent, measured along the envelope. This is the one that opens
	// the corner: with tip_out = 0 the far end sits ON the corner, the Bezier
	// collapses onto the leg it was meant to replace, and tip_in alone changes
	// the response by ~1e-15
	double tip_out = 0.0;

public:
	static const int NumParameters = 9;

	// whether this shape cuts the corner at zero stress at all, by either
	// mechanism
	inline bool hasChamfer() const {
		return chamfer > 0.0 || chamfer_rs > 0.0 || chamfer_dpe > 0.0;
	}
	// whether it cuts the reload/envelope junction. tip_out alone is enough: a
	// cut that only runs PAST the target still delays the engagement
	inline bool hasTip() const { return tip_in > 0.0 || tip_out > 0.0; }

	// beta selects the shape of the loop (0 = flag-shaped, 0.5 = peak-oriented,
	// 1 = classical plasticity), gamma how much permanent strain the closing
	// crack gives back. Both in [0, 1]
	static ASDHysteretic1DShape fromBetaGamma(
		double beta, double gamma, double chamfer_bias, double beta_bias,
		double smooth, double tip_in, double tip_out);
	// the same beta, with the corner cut given as a fraction of EACH leg.
	// Needed on the plastic half of the dial, where gammaScale() closes the
	// gamma gate and the corner cannot be rounded from there
	static ASDHysteretic1DShape fromBetaCut(
		double beta, double cut_in, double cut_out, double beta_bias,
		double smooth, double tip_in, double tip_out);
	// the same beta, with the corner cut given by its two ends
	static ASDHysteretic1DShape fromBetaRsDpe(
		double beta, double rs, double dpe, double beta_bias,
		double smooth, double tip_in, double tip_out);
	// how much of the requested gamma survives at this beta. The crack closing
	// back is a PINCHING mechanism, so it belongs to the half of the dial
	// between the flag-shaped and the peak-oriented loop and to nowhere else:
	// 4u(1-u) with u = min(2 beta, 1), identically zero from the
	// peak-oriented law onwards
	static double gammaScale(double beta);
	// slide the reloading vertex off the beta diagonal at fixed height: a3 and
	// a4 already ARE the pinchX/pinchY of the classical models, and all
	// fromBetaGamma does is derive one from the other
	static double betaBiasedA3(double a3_diagonal, double beta_bias);
	// the two cut fractions (incoming, outgoing) of a chamfer of size g. Built
	// so that neither can leave [0, 1], which is why no clamp is needed
	// anywhere: the one being lengthened goes from g toward 1 (never past its
	// neighbouring breakpoint) and the one being shortened from g toward 0
	static void chamferLegs(double g, double bias, double& in, double& out);

	// bounds check. On failure returns false and points 'bad' at the name of
	// the offending field and 'lo'/'hi' at its admissible range
	bool validate(const char** bad, double& lo, double& hi) const;

	// serialization
	void serialize(Vector& data, int& pos) const;
	void deserialize(const Vector& data, int& pos);
};

/**
A branch: a polyline sorted by strain, plus the strain at which it reaches zero
stress (its own first leg extended if it does not cross zero within its span).

FIXED CAPACITY, AND THE BOUND IS DERIVED. A branch is one birth point, a rounded
vertex and a rounded tip - SmoothSegments + 1 points each - and, once the closure
corner is cut, the points of the branch that survive the cut plus the rounded
chamfer. That is 40 at SmoothSegments = 12, measured on the calibrated setting.
The capacity is larger, and add() refuses rather than overflowing, so a change
that outgrows the bound fails loudly instead of corrupting memory.
*/
class ASDHysteretic1DBranch
{
public:
	// How many sub-segments a rounded corner is tessellated into. The curve is
	// approximated rather than carried symbolically, so that every piece of
	// machinery downstream - evaluation, the zero crossing, the free energy,
	// the IMPL-EX freeze - keeps working on a polyline exactly as before. The
	// error of a quadratic Bezier against its own chords falls as 1/n^2; at 12
	// it is under 1e-3 of the corner's own size, far below anything else the
	// model resolves.
	static const int SmoothSegments = 12;
	static const int MaxPoints = 4 * (SmoothSegments + 1) + 4;

public:
	double x[MaxPoints];
	double y[MaxPoints];
	int n = 0;
	double c = 0.0;

public:
	inline void clear() { n = 0; c = 0.0; }
	// false if the capacity is exhausted, which the bound above says cannot
	// happen
	inline bool add(double px, double py) {
		if (n >= MaxPoints)
			return false;
		x[n] = px;
		y[n] = py;
		++n;
		return true;
	}
	// (sigma, k) at eps, on the leg that holds eps_leg. eps_leg is the strain
	// that SELECTS the leg: the trial strain in the implicit solution, the
	// committed one under IMPL-EX - that is how the leg switch gets frozen,
	// like every other on/off state. Outside the branch's span the nearest leg
	// is extrapolated.
	void eval(double eps, double eps_leg, double& sigma, double& k) const;
};

class ASDHysteretic1DMaterial : public UniaxialMaterial, public IMPLEXObject
{
public:
	// sub-classes
	using HardeningLawPoint = ASDHardeningLawPoint;
	using HardeningLawType = ASDHardeningLawType;
	using HardeningLawPointComponent = ASDHardeningLawPointComponent;
	using HardeningLaw = ASDHardeningLaw;
	using BranchType = ASDHysteretic1DBranchType;
	using Shape = ASDHysteretic1DShape;
	using Branch = ASDHysteretic1DBranch;

	// Everything that places the reloading vertex and the zero crossing, with
	// no polyline and no tessellation, so that the engagement threshold can
	// reach the same numbers as the branch itself BY CONSTRUCTION instead of by
	// rebuilding the branch and reading its last point. That is how the two
	// used to drift apart.
	struct ReloadGeometry {
		const Shape* sh = 0;
		double tx = 0.0;
		double ty = 0.0;
		double vx = 0.0;
		double vy = 0.0;
		double c = 0.0;
	};

	// Everything the implicit pass of compute() can write. Saved and restored
	// around the non-destructive measurement of the IMPL-EX error, so that a
	// step that is measured and then rejected is a step that never happened.
	// It has to be COMPLETE: revertToLastCommit() does not restore side_commit,
	// which is a committed quantity the implicit pass overwrites.
	struct TrialState {
		double xt = 0.0;
		double xc = 0.0;
		BranchType bt_p = BranchType::ReloadCross;
		double pp_x = 0.0;
		double pp_y = 0.0;
		BranchType bt_n = BranchType::ReloadCross;
		double pn_x = 0.0;
		double pn_y = 0.0;
		double xb_p = 0.0;
		double xb_n = 0.0;
		int side_commit = 1;
		double strain = 0.0;
		double stress = 0.0;
		double stress_eff = 0.0;
		double C = 0.0;
		double dt_bar = 0.0;
		double dc_bar = 0.0;
		double anchor = 0.0;
	};

public:
	// life-cycle
	ASDHysteretic1DMaterial(
		int _tag,
		double _E,
		double _eta,
		bool _implex,
		bool _implex_control,
		bool _implex_abort_on_error,
		double _implex_error_tolerance,
		double _implex_time_reduction_limit,
		double _implex_alpha,
		bool _tangent,
		bool _auto_regularize,
		double _lch_ref,
		const HardeningLaw& _ht,
		const HardeningLaw& _hc,
		const Shape& _shape_t,
		const Shape& _shape_c,
		const std::vector<double>& _ls_t = std::vector<double>(),
		const std::vector<double>& _ls_c = std::vector<double>());
	ASDHysteretic1DMaterial();
	~ASDHysteretic1DMaterial();

	// info
	const char* getClassType(void) const { return "ASDHysteretic1DMaterial"; }

	// set strain
	int setTrialStrain(double v, double r = 0.0);

	// get state
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);

	// handle state
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// copy and others...
	UniaxialMaterial* getCopy(void);
	void Print(OPS_Stream& s, int flag = 0);

	// send/recv self
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);
	double getEnergy(void);

	// IMPL-EX error control
	double computeImplexErrorMetric(void);
	double implexTimeRatio(void) const;

private:
	// internal computation
	int compute(bool do_implex, bool do_tangent);
	// the IMPL-EX error metric, |s_implex - s_implicit| / s_ref, with s_ref the
	// largest stress this material can carry. The same quantity every material
	// of the family measures - see ASDConcrete1DMaterial for why it is not a
	// difference of damage (and here there is no damage to difference anyway)
	double implexStressGap(double delivered, double stress_implicit) const;
	// save/restore everything the implicit pass writes
	void saveTrialState(TrialState& x) const;
	void restoreTrialState(const TrialState& x);

	// --- branch geometry, mirroring the reference implementation ---------- //

	// (tx, ty, p) of one side: the envelope point at the side's current
	// measure (the first backbone corner while virgin) and the strain at which
	// the elastic line through it reaches zero stress, i.e. the side's FULL
	// plastic strain. All three SIGNED with the side
	void target(int side, double xt, double xc,
		double& tx, double& ty, double& p) const;
	// the shape a reloading branch of 'side' uses: the one of the side the
	// excursion STARTED from - the opposite one for a branch born at a zero
	// crossing, this one for a branch born at a reversal
	const Shape& reloadShape(int side, BranchType bt) const;
	// geometry of the current branch of 'side'. k_in is the slope the OTHER
	// side's branch has at the zero crossing: the reference of a negative a3
	// and the cap that forbids a stiffening kink at the closure point. Only
	// meaningful for ReloadCross; have_k_in = false disables it
	void branch(int side, BranchType bt, double px, double py,
		double xt, double xc, bool have_k_in, double k_in, Branch& out) const;
	// (sx, sy, c) of an UNLOADING branch: its start and its closure
	void unloadLeg(int side, BranchType bt, double px, double py,
		double xt, double xc, double& sx, double& sy, double& c) const;
	// its slope, without building the polyline. It has a single leg, so no
	// evaluation point is needed
	double unloadSlope(int side, BranchType bt, double px, double py,
		double xt, double xc) const;
	// the arithmetic of a RELOADING branch: no polyline, no tessellation
	ReloadGeometry reloadGeometry(int side, BranchType bt, double px, double py,
		double xt, double xc, bool have_k_in, double k_in) const;
	// the points that replace the bare target when the junction is cut.
	// Returns how many were written (0 = no cut)
	int tipPoints(int side, const HardeningLaw& law, const Shape& sh,
		double vx, double vy, double tx, double ty, double xtol,
		double* ox, double* oy) const;
	// the leg V -> T the junction cut measures itself against. The SINGLE
	// place that decides whether a reloading branch has a tip, so that the
	// branch and the engagement threshold can never disagree about it
	bool tipSpan(const HardeningLaw& law, const Shape& sh,
		double vx, double tx, double xtol, double& span) const;
	// strain at which a reloading branch stops running below the envelope: the
	// target when there is no tip, the far end of the tip when there is.
	// UNSIGNED, in the units of the side's own measure
	double reloadEnd(int side, BranchType bt, double px, double py,
		double xt, double xc) const;
	// the point next to the branch's END that sits at cx, which is the
	// direction the leg comes from or goes to
	bool neighbour(const Branch& br, double cx, double xtol,
		double& nx, double& ny) const;
	// (g_in, g_out) from the two DIRECT chamfer parameters
	void chamferDirect(const Shape& sh, double cx, double p_out_x, double xtol,
		double& g_in, double& g_out) const;
	// cut the corner both branches make at (cx, 0): br_in closes there,
	// br_out is born there. Both get the same A -> B segment and the same new
	// zero crossing, so the response stays continuous across the side change,
	// which happens on it
	void chamfer(const Branch& br_in, const Branch& br_out, double cx,
		double g, double xtol, double bias, const Shape* sh,
		Branch& new_in, Branch& new_out) const;
	// (corner, branch) of an unloading branch as it is actually walked: the
	// corner is where the branch reaches zero stress ON ITS OWN, which is
	// where the opposite side's reloading branch is born; the branch that
	// comes back has the corner cut, so its own c is where the stress actually
	// changes sign
	void unloadPath(int side, BranchType bt, double px, double py,
		double xt, double xc, double& corner, Branch& cut) const;
	// |int sigma deps| along a branch, from its zero-stress end to eps. Exact:
	// the branch is a polyline
	double area(const Branch& br, double eps) const;

	// --- responses -------------------------------------------------------- //
	Vector getHardeningLawVector(HardeningLawType ltype, HardeningLawPointComponent c) const;
	const Vector& getStrainMeasure() const;
	const Vector& getEquivalentPlasticStrain() const;
	const Vector& getClosureStrain() const;
	const Vector& getStiffnessReduction() const;
	const Vector& getImplexError() const;
	const Vector& getImplexStress() const;
	const Vector& getLimitStateRatio() const;
	const Vector& getTimeIncrements() const;
	const Vector& getShapeVector(int side) const;
	const Vector& getBranchState() const;
	// the energy recoverable by unloading to zero, integrated along the path
	// the material would ACTUALLY walk down - chamfer included
	const Vector& getFreeEnergy() const;

private:
	// Young's modulus
	double E = 0.0;
	// Viscosity for the rate-dependent update of the strain measures
	double eta = 0.0;
	// True = use the IMPL-EX algorithm
	bool implex = false;
	// True = measure the IMPL-EX error at every setTrialStrain and publish it.
	// The LEGACY path: it pays one implicit solve per Newton iteration, at
	// strains that are not equilibrated yet. The measurement the error control
	// runs on is the one CTestImplexWrapper asks for, once per step, through
	// computeImplexErrorMetric()
	bool implex_control = false;
	// True = let this material fail the step on its own when the error is over
	// tolerance. OFF by default and it should stay off: a material that aborts
	// owns a policy it cannot see the analysis to choose
	bool implex_abort_on_error = false;
	// Maximum allowed IMPL-EX error (default = 5%)
	double implex_error_tolerance = 0.05;
	// Minimum allowed time step reduction factor under which the legacy
	// in-material gate stops firing
	double implex_time_redution_limit = 0.01;
	// Scale factor for the implex extrapolation
	double implex_alpha = 1.0;
	// True = numerical tangent, False (default) = the active leg's slope,
	// which is the exact algorithmic operator here
	bool tangent = false;
	// True = regularize the fracture energy with the element's characteristic
	// length
	bool auto_regularize = true;
	bool regularization_done = false;
	double lch = 1.0;
	double lch_ref = 1.0;
	// the two backbones
	HardeningLaw ht;
	HardeningLaw hc;
	// the shape of the unload/reload path of each side (all zeros = the plain
	// bilinear peak-oriented law)
	Shape shape_t;
	Shape shape_c;
	// USER-DEFINED LIMIT STATES, one ascending list per side, in the material's
	// own STRAIN units - not as multiples of a yield strain. Deliberately not
	// normalized in here: "yield" has no single definition on this backbone.
	// strainAtOnsetOfCrack() is the start of SOFTENING, not the yield point, and
	// points()[1] is the yield point only if the user gave a bilinear start - on
	// a finely sampled curve it is a tiny elastic sub-step and the normalizer
	// would be wrong by an order of magnitude, silently. Whoever writes the input
	// knows the yield value and can multiply there. Empty = the feature is off
	std::vector<double> ls_t;
	std::vector<double> ls_c;
	// state variables - the equivalent strain measures (monotone, the only
	// extrapolated quantities)
	double xt = 0.0;
	double xt_commit = 0.0;
	double xt_commit_old = 0.0;
	double xc = 0.0;
	double xc_commit = 0.0;
	double xc_commit_old = 0.0;
	// on/off states, NEVER extrapolated and frozen at commit under IMPL-EX:
	// the branch type of each side and the point it was born at. The virgin
	// state is a reload from the origin to the first corner, i.e. the elastic
	// line
	BranchType bt_p = BranchType::ReloadCross;
	double pp_x = 0.0;
	double pp_y = 0.0;
	BranchType bt_n = BranchType::ReloadCross;
	double pn_x = 0.0;
	double pn_y = 0.0;
	BranchType bt_p_commit = BranchType::ReloadCross;
	double pp_x_commit = 0.0;
	double pp_y_commit = 0.0;
	BranchType bt_n_commit = BranchType::ReloadCross;
	double pn_x_commit = 0.0;
	double pn_y_commit = 0.0;
	// the side's own measure AT THE BIRTH of its reloading branch. Another
	// frozen switch, not an extrapolated quantity: it is what the reload AIMS
	// AT, and freezing it is what lets the strain run past the target over a
	// tip cut while the measure keeps tracking the strain. Without a tip it
	// always equals the live measure, since a side's measure only grows when
	// that side engages the envelope, which ends the reload
	double xb_p = 0.0;
	double xb_n = 0.0;
	double xb_p_commit = 0.0;
	double xb_n_commit = 0.0;
	// the active branch of the last implicit solution
	int side_commit = 1;
	// state variables - implex time bookkeeping
	double dtime_n = 0.0;
	double dtime_n_commit = 0.0;
	double dtime_0 = 0.0;
	bool dtime_is_user_defined = false;
	bool commit_done = false;
	double implex_error = 0.0;
	// strain, stress and tangent
	double strain = 0.0;
	double strain_commit = 0.0;
	double stress = 0.0;
	double stress_commit = 0.0;
	// WHAT THE STEP DELIVERED TO THE ELEMENT, kept because it is otherwise
	// unobservable: under IMPL-EX commitState re-solves implicitly and installs
	// that answer over `stress`, so from outside - a recorder, a test - the
	// extrapolated stress is gone by the time anyone can ask. It is what the
	// error metric measures the distance from. On the implicit path there is no
	// extrapolation and it equals `stress`
	double stress_implex = 0.0;
	// = stress: with no damage the nominal and the effective stress coincide.
	// Kept so that the output surface matches the rest of the family
	double stress_eff = 0.0;
	double stress_eff_commit = 0.0;
	double C = 0.0;
	// output only: 1 - k/E of the current leg of each side. NOT a damage
	// variable - nothing is stored and nothing is irreversible - it is the
	// instantaneous stiffness reduction of the leg the side is on
	double dt_bar = 0.0;
	double dc_bar = 0.0;
	// output only: the zero-stress crossing of the active branch
	double anchor = 0.0;
	double energy = 0.0;
};

#endif // ASDHysteretic1DMaterial_h
