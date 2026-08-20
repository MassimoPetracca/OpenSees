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
// See ASDHysteretic1DMaterial.h for the model.

#include <ASDHysteretic1DMaterial.h>
#include <Channel.h>
#include <ID.h>
#include <OPS_Globals.h>
#include <Information.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <Element.h>
#include <MaterialResponse.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

// anonymous namespace for utilities
namespace {

	enum ErrorCodes {
		EC_Generic = -1,
		EC_IMPLEX_Error_Control = -10
	};

	// The two exact ends of the beta dial have a leg of ZERO stiffness - a
	// plateau at zero stress for the flag-shaped loop, a plateau at the target
	// stress for the perfectly plastic one. They are geometrically right and
	// numerically unusable, so the dial is squeezed into [m, 1-m].
	const double SHAPE_MARGIN = 0.05;

	// How much of the corner gamma = 1 is allowed to cut. The geometric limit
	// is 1 (the chamfer then replaces both legs entirely), but a large cut
	// makes the unloading a SECANT far softer than E, and the dissipation on a
	// reloading leg, sigma*(1 - k/E_eff)*deps, turns negative as soon as that
	// leg is steeper than the secant. So this cap is a THERMODYNAMIC one and it
	// was measured, not guessed: at the binding case - the flag end, where the
	// reloading arrives on the target at the full E - a cut of 0.50 already
	// leaves a negative-work fraction that CONVERGES under refinement, i.e. a
	// real violation however small, while at 0.45 every cell of the
	// (beta, gamma) square decreases under refinement, which is the signature
	// of a quadrature artifact and of nothing else.
	const double CHAMFER_MAX = 0.45;

	/**
	Converts a string into a vector of doubles using whitespace as delimiter
	*/
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

	/**
	global parameters storage (shared convention with the other ASD materials)
	*/
	class GlobalParameters {
	private:
		double max_error = 0.0;
		double avg_error = 0.0;
		int avg_counter = 0;
	private:
		GlobalParameters() = default;
		GlobalParameters(const GlobalParameters&) = delete;
		GlobalParameters& operator = (const GlobalParameters&) = delete;
	public:
		static GlobalParameters& instance() {
			static GlobalParameters _instance;
			return _instance;
		}
		inline double getMaxError() const { return max_error; }
		inline void setMaxError(double x) { max_error = x; }
		inline double getAverageError() {
			if (avg_counter > 0) {
				avg_error /= static_cast<double>(avg_counter);
				avg_counter = 0;
			}
			return avg_error;
		}
		inline void accumulateAverageError(double x) {
			avg_error += x;
			++avg_counter;
		}
		inline void setAverageError(double x) {
			avg_error = x;
			avg_counter = 0;
		}
	};

	inline double clamp01(double x) {
		return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
	}

	/**
	The candidate points of a branch, before they are sorted and de-duplicated
	into it. A plain fixed-capacity buffer for the same reason the branch is
	one: this is the hottest allocation of the model, one per branch per Gauss
	point per iteration, and it does not need to exist.
	*/
	struct PointBuffer {
		static const int Capacity = ASDHysteretic1DBranch::MaxPoints;
		double x[Capacity];
		double y[Capacity];
		int n = 0;
		inline void clear() { n = 0; }
		inline void add(double px, double py) {
			if (n < Capacity) {
				x[n] = px;
				y[n] = py;
				++n;
			}
		}
		inline void add(const double* px, const double* py, int count) {
			for (int i = 0; i < count; ++i)
				add(px[i], py[i]);
		}
	};

	/**
	Sort the candidates by strain, drop the degenerate ones, and store them in
	the branch together with its zero crossing.

	THE SORT MUST BE STABLE, and it is not a safety net that could be dropped.
	The points of a branch do arrive monotone, but the chamfer hands over the
	KEPT points of a branch followed by the cut that replaces its corner, and
	for the branch BORN at the corner those two come in the wrong order. Where
	two points share a strain the FIRST one given wins, which is why an
	insertion sort and not std::sort.
	*/
	void polyline(const PointBuffer& in, double c, double xtol, double E,
		ASDHysteretic1DBranch& out)
	{
		// stable insertion sort into the output arrays
		double sx[PointBuffer::Capacity];
		double sy[PointBuffer::Capacity];
		int m = 0;
		for (int i = 0; i < in.n; ++i) {
			int j = m;
			while (j > 0 && sx[j - 1] > in.x[i]) {
				sx[j] = sx[j - 1];
				sy[j] = sy[j - 1];
				--j;
			}
			sx[j] = in.x[i];
			sy[j] = in.y[i];
			++m;
		}
		// de-duplicate
		out.clear();
		out.c = c;
		if (m > 0) {
			out.add(sx[0], sy[0]);
			for (int i = 1; i < m; ++i) {
				if (sx[i] - out.x[out.n - 1] > xtol)
					out.add(sx[i], sy[i]);
			}
		}
		if (out.n < 2) {
			// fully degenerate branch: fall back on the elastic line through
			// (c, 0). The span is irrelevant, evaluation extrapolates anyway
			out.clear();
			out.c = c;
			out.add(c, 0.0);
			out.add(c + 1.0, E);
		}
	}

	/**
	Round the corner instead of cutting it straight, and say where it crosses
	zero stress.

	The cut A -> B is replaced by the QUADRATIC BEZIER on A, a control point C
	and B, with C = (1 - smooth)*midpoint(A, B) + smooth*corner, so smooth = 0
	puts the control point on the chord and the curve IS the straight cut, while
	smooth = 1 puts it on the corner and the curve leaves A tangent to the
	incoming leg and reaches B tangent to the outgoing one - C1 with both, where
	the chamfer has two kinks.

	A quadratic Bezier and not a NURBS on purpose. The whole admissibility
	argument of the chamfer is that the cut lies between the two legs it joins,
	which here comes for free: the curve stays in the triangle A C B and its
	tangent, 2[(1-t)(C-A) + t(B-C)], is a convex combination of the two leg
	directions, so its slope can never leave the interval they span.

	The corner sits AT zero stress, so the stress along the curve is a plain
	quadratic and its root is closed-form - no iteration, and the new closure
	point stays exact.
	*/
	int roundCorner(double ax, double ay, double cx, double bx, double by,
		double smooth, double* ox, double* oy, double& c)
	{
		const int n = ASDHysteretic1DBranch::SmoothSegments;
		double s = clamp01(smooth);
		double ccx = (1.0 - s) * 0.5 * (ax + bx) + s * cx;
		double ccy = (1.0 - s) * 0.5 * (ay + by);   // the corner's own stress is 0
		for (int i = 0; i <= n; ++i) {
			double t = static_cast<double>(i) / static_cast<double>(n);
			double u = 1.0 - t;
			ox[i] = u * u * ax + 2.0 * t * u * ccx + t * t * bx;
			oy[i] = u * u * ay + 2.0 * t * u * ccy + t * t * by;
		}
		// sigma(t) = 0, closed form; fall back on the chord if it degenerates
		double qa = ay - 2.0 * ccy + by;
		double qb = 2.0 * (ccy - ay);
		double qc = ay;
		bool have_t0 = false;
		double t0 = 0.0;
		if (std::abs(qa) < 1.0e-300) {
			if (qb != 0.0) {
				t0 = -qc / qb;
				have_t0 = true;
			}
		}
		else {
			double disc = qb * qb - 4.0 * qa * qc;
			if (disc >= 0.0) {
				double r = std::sqrt(disc);
				double roots[2] = { (-qb + r) / (2.0 * qa), (-qb - r) / (2.0 * qa) };
				for (int k = 0; k < 2; ++k) {
					if (-1.0e-12 <= roots[k] && roots[k] <= 1.0 + 1.0e-12) {
						t0 = clamp01(roots[k]);
						have_t0 = true;
						break;
					}
				}
			}
		}
		if (!have_t0) {
			c = (by != ay) ? (ax - ay * (bx - ax) / (by - ay)) : cx;
		}
		else {
			double u = 1.0 - t0;
			c = u * u * ax + 2.0 * t0 * u * ccx + t0 * t0 * bx;
		}
		return n + 1;
	}

	/**
	Round (or cut straight) the corner where the reload meets the envelope.

	The third corner of the excursion, and the simplest of the three: the two
	ends are already placed by tip_in and tip_out, so unlike the vertex there is
	no extent to invent, and unlike the closure corner there is no zero crossing
	to solve for - this one sits ON the envelope, nowhere near sigma = 0.

	The tangent is again a convex combination of the two legs, so the rounded
	path stays inside the corner: the response can only fall BELOW the envelope
	here, never above it, whatever smooth is.
	*/
	int roundTip(double ax, double ay, double cx, double cy, double bx, double by,
		double smooth, double* ox, double* oy)
	{
		double s = clamp01(smooth);
		if (s <= 0.0) {
			ox[0] = ax; oy[0] = ay;
			ox[1] = bx; oy[1] = by;
			return 2;
		}
		const int n = ASDHysteretic1DBranch::SmoothSegments;
		double ccx = (1.0 - s) * 0.5 * (ax + bx) + s * cx;
		double ccy = (1.0 - s) * 0.5 * (ay + by) + s * cy;
		for (int i = 0; i <= n; ++i) {
			double t = static_cast<double>(i) / static_cast<double>(n);
			double u = 1.0 - t;
			ox[i] = u * u * ax + 2.0 * t * u * ccx + t * t * bx;
			oy[i] = u * u * ay + 2.0 * t * u * ccy + t * t * by;
		}
		return n + 1;
	}

	/**
	Round the RELOADING vertex, the other corner of the excursion.

	There is no chamfer here to blend with, so smooth sets the EXTENT instead of
	the blend: half of each adjacent leg at smooth = 1, nothing at 0, and the
	control point is the vertex itself. Half is the natural maximum - it is as
	far as both cuts can go without either reaching a neighbour - so the
	guarantee that a rounding never overruns a leg holds with no clamp.

	Writes just the vertex when there is nothing to round, so the caller is
	unchanged.
	*/
	int roundVertex(double sx, double sy, double vx, double vy,
		double tx, double ty, double smooth, double xtol,
		double* ox, double* oy)
	{
		double s = clamp01(smooth);
		if (s <= 0.0 || std::abs(vx - sx) <= xtol || std::abs(tx - vx) <= xtol) {
			ox[0] = vx; oy[0] = vy;
			return 1;
		}
		const int n = ASDHysteretic1DBranch::SmoothSegments;
		double h = 0.5 * s;
		double a0 = vx + h * (sx - vx), a1 = vy + h * (sy - vy);
		double b0 = vx + h * (tx - vx), b1 = vy + h * (ty - vy);
		for (int i = 0; i <= n; ++i) {
			double t = static_cast<double>(i) / static_cast<double>(n);
			double u = 1.0 - t;
			ox[i] = u * u * a0 + 2.0 * t * u * vx + t * t * b0;
			oy[i] = u * u * a1 + 2.0 * t * u * vy + t * t * b1;
		}
		return n + 1;
	}

	// the largest number of points any of the three roundings writes
	const int MaxRoundPoints = ASDHysteretic1DBranch::SmoothSegments + 1;
}

// ---------------------------------------------------------------------------
// ASDHysteretic1DShape
// ---------------------------------------------------------------------------

double ASDHysteretic1DShape::gammaScale(double beta)
{
	double u = 2.0 * beta;
	if (u > 1.0) u = 1.0;
	return 4.0 * u * (1.0 - u);
}

double ASDHysteretic1DShape::betaBiasedA3(double a3_diagonal, double beta_bias)
{
	double t = beta_bias;
	if (t > 1.0) t = 1.0;
	if (t < -1.0) t = -1.0;
	if (t == 0.0)
		return a3_diagonal;
	return a3_diagonal + t * ((t > 0.0) ? (1.0 - a3_diagonal) : (a3_diagonal + 1.0));
}

void ASDHysteretic1DShape::chamferLegs(double g, double bias, double& in, double& out)
{
	double s = std::abs(bias);
	if (s > 1.0) s = 1.0;
	double long_ = g + (1.0 - g) * s;
	double short_ = g * (1.0 - s);
	if (bias >= 0.0) {
		in = long_;
		out = short_;
	}
	else {
		in = short_;
		out = long_;
	}
}

ASDHysteretic1DShape ASDHysteretic1DShape::fromBetaGamma(
	double beta, double gamma, double chamfer_bias, double beta_bias,
	double smooth, double tip_in, double tip_out)
{
	// The vertex is placed at the fraction b of the segment joining the two
	// extreme useful positions - the target's own elastic crossing (a flag) and
	// the point where the incoming slope carried past zero stress reaches the
	// target's stress (classical plasticity). The MIDPOINT of that segment lies
	// exactly on the peak-oriented chord, hence b = 1/2 is the previous law bit
	// for bit. Since the vertex sits at the height b of the span, a4 = b, and
	// the piecewise a3 that puts it at the fraction b of the FULL span is
	// (1 - 2b)/(1 - b) below the chord and -(2b - 1)/b above it.
	double m = SHAPE_MARGIN;
	double b = m + beta * (1.0 - 2.0 * m);
	double a3 = (b <= 0.5) ? ((1.0 - 2.0 * b) / (1.0 - b)) : (-(2.0 * b - 1.0) / b);
	ASDHysteretic1DShape sh;
	sh.a3 = betaBiasedA3(a3, beta_bias);
	sh.a4 = b;
	sh.chamfer = CHAMFER_MAX * gamma * gammaScale(beta);
	sh.chamfer_bias = chamfer_bias;
	sh.smooth = smooth;
	sh.tip_in = tip_in;
	sh.tip_out = tip_out;
	return sh;
}

ASDHysteretic1DShape ASDHysteretic1DShape::fromBetaCut(
	double beta, double cut_in, double cut_out, double beta_bias,
	double smooth, double tip_in, double tip_out)
{
	// The corner cut given as a fraction of EACH leg, one per leg. Stored as
	// the equivalent (chamfer, chamfer_bias), so this is not a third
	// mechanism: it is a naming of the one chamferLegs() already implements.
	// Inverting it, with lo/hi the smaller and larger of the two:
	//     chamfer = lo / (lo + 1 - hi),   bias = +-(1 - lo/chamfer)
	// the sign being positive when the incoming leg is the longer cut.
	double a = clamp01(cut_in);
	double b = clamp01(cut_out);
	ASDHysteretic1DShape base = fromBetaGamma(beta, 0.0, 0.0, beta_bias, 0.0, 0.0, 0.0);
	ASDHysteretic1DShape sh;
	sh.a3 = base.a3;
	sh.a4 = base.a4;
	sh.tip_in = tip_in;
	sh.tip_out = tip_out;
	if (a <= 0.0 || b <= 0.0) {
		// a cut that reaches zero on one leg is a cut along the other one,
		// i.e. no corner removed at all
		return sh;
	}
	double lo = (a >= b) ? b : a;
	double hi = (a >= b) ? a : b;
	double den = lo + 1.0 - hi;
	double g = (den <= 0.0) ? 1.0 : (lo / den);
	double s = (g <= 0.0) ? 0.0 : (1.0 - lo / g);
	sh.chamfer = g;
	sh.chamfer_bias = (a >= b) ? s : -s;
	sh.smooth = smooth;
	return sh;
}

ASDHysteretic1DShape ASDHysteretic1DShape::fromBetaRsDpe(
	double beta, double rs, double dpe, double beta_bias,
	double smooth, double tip_in, double tip_out)
{
	ASDHysteretic1DShape base = fromBetaGamma(beta, 0.0, 0.0, beta_bias, 0.0, 0.0, 0.0);
	ASDHysteretic1DShape sh;
	sh.a3 = base.a3;
	sh.a4 = base.a4;
	sh.chamfer_rs = rs;
	sh.chamfer_dpe = dpe;
	sh.smooth = smooth;
	sh.tip_in = tip_in;
	sh.tip_out = tip_out;
	return sh;
}

bool ASDHysteretic1DShape::validate(const char** bad, double& lo, double& hi) const
{
	// the shape parameters are dimensionless fractions
	struct Item { const char* name; double value; double lo; double hi; };
	const Item items[NumParameters] = {
		{ "a3", a3, -1.0, 1.0 },
		{ "a4", a4, 0.0, 1.0 },
		{ "chamfer", chamfer, 0.0, 1.0 },
		{ "chamferBias", chamfer_bias, -1.0, 1.0 },
		{ "chamferRs", chamfer_rs, 0.0, 1.0 },
		{ "chamferDpe", chamfer_dpe, 0.0, 1.0 },
		{ "smooth", smooth, 0.0, 1.0 },
		{ "tipIn", tip_in, 0.0, 1.0 },
		{ "tipOut", tip_out, 0.0, 1.0 }
	};
	for (int i = 0; i < NumParameters; ++i) {
		if (!(items[i].lo <= items[i].value && items[i].value <= items[i].hi)) {
			*bad = items[i].name;
			lo = items[i].lo;
			hi = items[i].hi;
			return false;
		}
	}
	return true;
}

void ASDHysteretic1DShape::serialize(Vector& data, int& pos) const
{
	data(pos++) = a3;
	data(pos++) = a4;
	data(pos++) = chamfer;
	data(pos++) = chamfer_bias;
	data(pos++) = chamfer_rs;
	data(pos++) = chamfer_dpe;
	data(pos++) = smooth;
	data(pos++) = tip_in;
	data(pos++) = tip_out;
}

void ASDHysteretic1DShape::deserialize(const Vector& data, int& pos)
{
	a3 = data(pos++);
	a4 = data(pos++);
	chamfer = data(pos++);
	chamfer_bias = data(pos++);
	chamfer_rs = data(pos++);
	chamfer_dpe = data(pos++);
	smooth = data(pos++);
	tip_in = data(pos++);
	tip_out = data(pos++);
}

// ---------------------------------------------------------------------------
// ASDHysteretic1DBranch
// ---------------------------------------------------------------------------

void ASDHysteretic1DBranch::eval(double eps, double eps_leg, double& sigma, double& k) const
{
	int i = 0;
	while (i < n - 2 && eps_leg > x[i + 1])
		++i;
	double xa = x[i], ya = y[i];
	double xb = x[i + 1], yb = y[i + 1];
	k = (yb - ya) / (xb - xa);
	sigma = ya + (eps - xa) * k;
}

// ---------------------------------------------------------------------------
// the command
// ---------------------------------------------------------------------------

void* OPS_ASDHysteretic1DMaterial(void)
{
	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 2) {
		opserr <<
			"uniaxialMaterial ASDHysteretic1D Error: Few arguments (< 2).\n"
			"uniaxialMaterial ASDHysteretic1D $tag $E "
			"-Te $Te -Ts $Ts -Ce $Ce -Cs $Cs "
			"<-shapeT $beta $gamma> <-shapeC $beta $gamma> "
			"<-shapeCutT $beta $cutIn $cutOut> <-shapeCutC ...> "
			"<-shapeRsDpeT $beta $rs $dpe> <-shapeRsDpeC ...> "
			"<-shapeRawT $a3 $a4 $chamfer $chamferBias $chamferRs $chamferDpe $smooth $tipIn $tipOut> <-shapeRawC ...> "
			"<-smoothT $s> <-smoothC $s> <-tipT $tipIn $tipOut> <-tipC $tipIn $tipOut> "
			"<-biasT $chamferBias $betaBias> <-biasC $chamferBias $betaBias> "
			"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-implexAbort> <-implexAlpha $alpha> "
			"<-eta $eta> <-tangent> <-autoRegularization $lch_ref> "
			"<-limitStates $x1 $x2 ...> <-limitStatesC $x1 $x2 ...>\n";
		return nullptr;
	}

	// numData
	int numData = 1;

	// data
	int tag;
	double E;
	bool implex = false;
	bool implex_control = false;
	bool implex_abort_on_error = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	double implex_alpha = 1.0;
	double eta = 0.0;
	bool tangent = false;
	bool auto_regularization = false;
	double lch_ref = 1.0;
	std::vector<double> Te, Ts, Td, Ce, Cs, Cd;
	// the user's limit states, ascending, in strain units. One list, applied to
	// both sides, unless the compressive one is given separately
	std::vector<double> LSt, LSc;
	bool has_LSc = false;

	// The shape is assembled from independent flags and only built at the end,
	// so their order does not matter: one of the four BASE flags writes
	// a3/a4/chamfer*, the modifiers overlay smooth, the tip and the two biases
	// whichever base was used.
	enum ShapeMode { SM_Raw = 0, SM_BetaGamma, SM_BetaCut, SM_BetaRsDpe };
	struct ShapeInput {
		ShapeMode mode = SM_Raw;
		double raw[ASDHysteretic1DShape::NumParameters] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double beta = 0.0;
		double p1 = 0.0;          // gamma, cutIn or rs
		double p2 = 0.0;          // -,     cutOut or dpe
		double chamfer_bias = 0.0;
		double beta_bias = 0.0;
		double smooth = 0.0;
		double tip_in = 0.0;
		double tip_out = 0.0;
		ASDHysteretic1DShape build() const {
			ASDHysteretic1DShape sh;
			switch (mode) {
			case SM_BetaGamma:
				sh = ASDHysteretic1DShape::fromBetaGamma(beta, p1, chamfer_bias, beta_bias, smooth, tip_in, tip_out);
				break;
			case SM_BetaCut:
				sh = ASDHysteretic1DShape::fromBetaCut(beta, p1, p2, beta_bias, smooth, tip_in, tip_out);
				break;
			case SM_BetaRsDpe:
				sh = ASDHysteretic1DShape::fromBetaRsDpe(beta, p1, p2, beta_bias, smooth, tip_in, tip_out);
				break;
			default:
				sh.a3 = raw[0];
				sh.a4 = raw[1];
				sh.chamfer = raw[2];
				sh.chamfer_bias = raw[3];
				sh.chamfer_rs = raw[4];
				sh.chamfer_dpe = raw[5];
				sh.smooth = raw[6];
				sh.tip_in = raw[7];
				sh.tip_out = raw[8];
				// the modifiers still win, so -smoothT after -shapeRawT is not
				// silently ignored
				if (smooth != 0.0) sh.smooth = smooth;
				if (tip_in != 0.0) sh.tip_in = tip_in;
				if (tip_out != 0.0) sh.tip_out = tip_out;
				if (chamfer_bias != 0.0) sh.chamfer_bias = chamfer_bias;
				break;
			}
			return sh;
		}
	};
	ShapeInput sit, sic;

	// get tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: invalid 'tag'.\n";
		return nullptr;
	}

	// get Elasticity arguments
	if (OPS_GetDouble(&numData, &E) != 0) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: invalid 'E'.\n";
		return nullptr;
	}
	if (E <= 0.0) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: invalid value for 'E' (" << E << "). It should be strictly positive.\n";
		return nullptr;
	}

	// utilities (code re-use)
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "uniaxialMaterial ASDHysteretic1D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "uniaxialMaterial ASDHysteretic1D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_list = [&numData](const char* variable, std::vector<double>& value) -> bool {
		// first try expanded list like {*}$the_list,
		// also used in python like *the_list
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
		// try Tcl list (it's a string after all...)
		if (value.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
			std::string list_string = OPS_GetString();
			if (!string_to_list_of_doubles(list_string, ' ', value)) {
				opserr << "uniaxialMaterial ASDHysteretic1D Error: cannot parse the '" << variable << "' list.\n";
				return false;
			}
		}
		return true;
	};
	// the n numbers of a shape flag, all mandatory
	auto lam_shape_numbers = [&numData](const char* variable, int count, double* out) -> bool {
		if (OPS_GetNumRemainingInputArgs() < count) {
			opserr << "uniaxialMaterial ASDHysteretic1D Error: '" << variable
				<< "' requires " << count << " values.\n";
			return false;
		}
		for (int i = 0; i < count; ++i) {
			if (OPS_GetDouble(&numData, &out[i]) < 0) {
				opserr << "uniaxialMaterial ASDHysteretic1D Error: failed to get value "
					<< (i + 1) << " of '" << variable << "'.\n";
				return false;
			}
		}
		return true;
	};

	// optional parameters
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "uniaxialMaterial ASDHysteretic1D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
				return nullptr;
			}
			if (!lam_optional_double("implexErrorTolerance", implex_error_tolerance))
				return nullptr;
			if (!lam_optional_double("implexTimeReductionLimit", implex_time_redution_limit))
				return nullptr;
		}
		else if (strcmp(value, "-implexAbort") == 0) {
			// the legacy behaviour: let the material fail the step by itself.
			// See ASDHysteretic1DMaterial.h on implex_abort_on_error
			implex_abort_on_error = true;
		}
		else if (strcmp(value, "-implexAlpha") == 0) {
			if (!lam_optional_double("alpha", implex_alpha))
				return nullptr;
		}
		else if (strcmp(value, "-eta") == 0) {
			if (!lam_optional_double("eta", eta))
				return nullptr;
		}
		else if (strcmp(value, "-tangent") == 0) {
			tangent = true;
		}
		else if (strcmp(value, "-autoRegularization") == 0) {
			auto_regularization = true;
			if (!lam_optional_double("lch_ref", lch_ref))
				return nullptr;
		}
		else if (strcmp(value, "-limitStates") == 0) {
			if (!lam_optional_list("limitStates", LSt))
				return nullptr;
		}
		else if (strcmp(value, "-limitStatesC") == 0) {
			if (!lam_optional_list("limitStatesC", LSc))
				return nullptr;
			has_LSc = true;
		}
		else if (strcmp(value, "-Te") == 0) {
			if (!lam_optional_list("Te", Te))
				return nullptr;
		}
		else if (strcmp(value, "-Ts") == 0) {
			if (!lam_optional_list("Ts", Ts))
				return nullptr;
		}
		else if (strcmp(value, "-Td") == 0) {
			if (!lam_optional_list("Td", Td))
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
		else if (strcmp(value, "-Cd") == 0) {
			if (!lam_optional_list("Cd", Cd))
				return nullptr;
		}
		// --- the shape of each side ---------------------------------------
		else if (strcmp(value, "-shapeT") == 0 || strcmp(value, "-shapeC") == 0) {
			ShapeInput& s = (value[6] == 'T') ? sit : sic;
			double v[2];
			if (!lam_shape_numbers(value, 2, v))
				return nullptr;
			s.mode = SM_BetaGamma;
			s.beta = v[0];
			s.p1 = v[1];
		}
		else if (strcmp(value, "-shapeCutT") == 0 || strcmp(value, "-shapeCutC") == 0) {
			ShapeInput& s = (value[9] == 'T') ? sit : sic;
			double v[3];
			if (!lam_shape_numbers(value, 3, v))
				return nullptr;
			s.mode = SM_BetaCut;
			s.beta = v[0];
			s.p1 = v[1];
			s.p2 = v[2];
		}
		else if (strcmp(value, "-shapeRsDpeT") == 0 || strcmp(value, "-shapeRsDpeC") == 0) {
			ShapeInput& s = (value[11] == 'T') ? sit : sic;
			double v[3];
			if (!lam_shape_numbers(value, 3, v))
				return nullptr;
			s.mode = SM_BetaRsDpe;
			s.beta = v[0];
			s.p1 = v[1];
			s.p2 = v[2];
		}
		else if (strcmp(value, "-shapeRawT") == 0 || strcmp(value, "-shapeRawC") == 0) {
			ShapeInput& s = (value[9] == 'T') ? sit : sic;
			s.mode = SM_Raw;
			if (!lam_shape_numbers(value, ASDHysteretic1DShape::NumParameters, s.raw))
				return nullptr;
		}
		else if (strcmp(value, "-smoothT") == 0 || strcmp(value, "-smoothC") == 0) {
			ShapeInput& s = (value[7] == 'T') ? sit : sic;
			if (!lam_optional_double(value, s.smooth))
				return nullptr;
		}
		else if (strcmp(value, "-tipT") == 0 || strcmp(value, "-tipC") == 0) {
			ShapeInput& s = (value[4] == 'T') ? sit : sic;
			double v[2];
			if (!lam_shape_numbers(value, 2, v))
				return nullptr;
			s.tip_in = v[0];
			s.tip_out = v[1];
		}
		else if (strcmp(value, "-biasT") == 0 || strcmp(value, "-biasC") == 0) {
			ShapeInput& s = (value[5] == 'T') ? sit : sic;
			double v[2];
			if (!lam_shape_numbers(value, 2, v))
				return nullptr;
			s.chamfer_bias = v[0];
			s.beta_bias = v[1];
		}
	}

	// check lists. There is no -fc preset here, on purpose: generating the
	// backbones from a strength also generates a DAMAGE, which this material
	// has to refuse (see below), so the preset of ASDConcrete1D cannot simply
	// be reused. Give the two backbones explicitly
	if (Te.size() < 1) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: 'Te' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Ts.size() != Te.size()) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Ts' (size = " <<
			static_cast<int>(Ts.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Ce.size() < 1) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: 'Ce' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Cs.size() != Ce.size()) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cs' (size = " <<
			static_cast<int>(Cs.size()) << ") lists should have the same size.\n";
		return nullptr;
	}

	// THE DAMAGE MUST BE ZERO, and it is refused rather than ignored: a
	// non-zero damage would silently reshape the backbone inside
	// ASDHardeningLaw::adjust(), so accepting it would change the envelope
	// while pretending nothing happened. This material has no damage variable -
	// the whole inelastic strain is plastic and the unilateral effect comes
	// from the branch topology
	for (int k = 0; k < 2; ++k) {
		const std::vector<double>& d = (k == 0) ? Td : Cd;
		const char* name = (k == 0) ? "Td" : "Cd";
		const std::vector<double>& e = (k == 0) ? Te : Ce;
		if (d.size() == 0)
			continue;
		if (d.size() != e.size()) {
			opserr << "uniaxialMaterial ASDHysteretic1D Error: '" << name << "' (size = " <<
				static_cast<int>(d.size()) << ") should have the same size as its strain list (size = " <<
				static_cast<int>(e.size()) << ").\n";
			return nullptr;
		}
		for (std::size_t i = 0; i < d.size(); ++i) {
			if (d[i] != 0.0) {
				opserr << "uniaxialMaterial ASDHysteretic1D Error: '" << name
					<< "' must be zero (or omitted): this material has no damage, the whole "
					"inelastic strain is plastic. A non-zero damage would also reshape the "
					"backbone through the hardening law's own adjustment.\n";
				return nullptr;
			}
		}
	}
	Td.assign(Te.size(), 0.0);
	Cd.assign(Ce.size(), 0.0);

	// build the hardening laws
	ASDHardeningLaw ht(tag, ASDHardeningLawType::Tension, E, Te, Ts, Td);
	if (!ht.isValid()) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: invalid tensile hardening law.\n";
		return nullptr;
	}
	ASDHardeningLaw hc(tag, ASDHardeningLawType::Compression, E, Ce, Cs, Cd);
	if (!hc.isValid()) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: invalid compressive hardening law.\n";
		return nullptr;
	}

	// build the shapes and check their bounds
	ASDHysteretic1DShape shape_t = sit.build();
	ASDHysteretic1DShape shape_c = sic.build();
	for (int k = 0; k < 2; ++k) {
		const ASDHysteretic1DShape& sh = (k == 0) ? shape_t : shape_c;
		const char* side = (k == 0) ? "shapeT" : "shapeC";
		const char* bad = "";
		double lo = 0.0, hi = 0.0;
		if (!sh.validate(&bad, lo, hi)) {
			opserr << "uniaxialMaterial ASDHysteretic1D Error: invalid value for '"
				<< side << "." << bad << "'. The shape parameters are dimensionless "
				"fractions, this one in [" << lo << ", " << hi << "].\n";
			return nullptr;
		}
	}

	// the limit states. One list serves both sides unless the compressive one was
	// given: the thresholds are ABSCISSAE, so an asymmetric backbone needs its own
	if (!has_LSc)
		LSc = LSt;
	for (int k = 0; k < 2; ++k) {
		const std::vector<double>& ls = (k == 0) ? LSt : LSc;
		const char* name = (k == 0) ? "limitStates" : "limitStatesC";
		for (std::size_t i = 0; i < ls.size(); ++i) {
			// STRICTLY ASCENDING AND POSITIVE, and refused rather than sorted: a
			// list out of order means the user meant something else, and two equal
			// thresholds would divide by zero in the interpolation between them
			if (!(ls[i] > 0.0)) {
				opserr << "uniaxialMaterial ASDHysteretic1D Error: '" << name
					<< "' must be positive, got " << ls[i] << " at position "
					<< static_cast<int>(i + 1) << ". They are ABSCISSAE of the "
					"strain measure, which is the peak |strain| reached on the side, "
					"so they are positive on both sides.\n";
				return nullptr;
			}
			if (i > 0 && !(ls[i] > ls[i - 1])) {
				opserr << "uniaxialMaterial ASDHysteretic1D Error: '" << name
					<< "' must be strictly ascending: " << ls[i - 1] << " then "
					<< ls[i] << ".\n";
				return nullptr;
			}
		}
	}

	// create the material
	UniaxialMaterial* instance = new ASDHysteretic1DMaterial(
		tag,
		E, eta,
		implex, implex_control, implex_abort_on_error,
		implex_error_tolerance, implex_time_redution_limit, implex_alpha,
		tangent,
		auto_regularization, lch_ref,
		ht, hc, shape_t, shape_c, LSt, LSc);
	if (instance == nullptr) {
		opserr << "uniaxialMaterial ASDHysteretic1D Error: failed to allocate a new material.\n";
		return nullptr;
	}
	return instance;
}

// ---------------------------------------------------------------------------
// life cycle
// ---------------------------------------------------------------------------

ASDHysteretic1DMaterial::ASDHysteretic1DMaterial(
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
	const std::vector<double>& _ls_t,
	const std::vector<double>& _ls_c)
	: UniaxialMaterial(_tag, MAT_TAG_ASDHysteretic1DMaterial)
	, E(_E)
	, eta(_eta)
	, implex(_implex)
	, implex_control(_implex_control)
	, implex_abort_on_error(_implex_abort_on_error)
	, implex_error_tolerance(_implex_error_tolerance)
	, implex_time_redution_limit(_implex_time_reduction_limit)
	, implex_alpha(_implex_alpha)
	, tangent(_tangent)
	, auto_regularize(_auto_regularize)
	, lch_ref(_lch_ref)
	, ht(_ht)
	, hc(_hc)
	, shape_t(_shape_t)
	, shape_c(_shape_c)
	, ls_t(_ls_t)
	, ls_c(_ls_c)
{
	C = E;
}

ASDHysteretic1DMaterial::ASDHysteretic1DMaterial()
	: UniaxialMaterial(0, MAT_TAG_ASDHysteretic1DMaterial)
{
}

ASDHysteretic1DMaterial::~ASDHysteretic1DMaterial()
{
}

// ---------------------------------------------------------------------------
// set strain
// ---------------------------------------------------------------------------

int ASDHysteretic1DMaterial::setTrialStrain(double v, double r)
{
	// return value
	int retval = 0;

	// this material point takes part in the current step, so it takes part in
	// the IMPL-EX error aggregate (see IMPLEXManager.h)
	implexTouch();

	// get characteristic length and perform regularization
	if (!regularization_done) {
		if (ops_TheActiveElement)
			lch = ops_TheActiveElement->getCharacteristicLength();
		regularization_done = true;
		if (auto_regularize) {
			ht.regularize(lch, lch_ref);
			hc.regularize(lch, lch_ref);
		}
	}

	// save dT
	if (!dtime_is_user_defined) {
		dtime_n = ops_Dt;
		if (!commit_done) {
			dtime_0 = dtime_n;
			dtime_n_commit = dtime_n;
		}
	}

	// if the user requested the numerical tangent and not the IMPL-EX (in
	// IMPL-EX the tangent coincides with the secant)...
	if (tangent && !implex) {
		double Cnum = 0.0;
		double PERT = (ht.strainTolerance() + hc.strainTolerance()) / 2.0;
		strain = v + PERT;
		retval = compute(true, false);
		if (retval < 0) return retval;
		Cnum = stress;
		strain = v;
		retval = compute(true, false);
		if (retval < 0) return retval;
		C = (Cnum - stress) / PERT;
	}
	else {
		strain = v;
		if (implex) {
			if (implex_control) {
				// LEGACY in-material measurement: an implicit solution first,
				// to have something to compare against, then the extrapolated
				// one, which is the answer that must reach the element. The
				// error control does not need this: CTestImplexWrapper
				// measures once per step, through computeImplexErrorMetric(),
				// instead of once per iteration
				int aux = side_commit;
				retval = compute(false, false);
				if (retval < 0) return retval;
				double stress_implicit = stress;
				// the implicit pass has overwritten the frozen switch: put it
				// back before extrapolating, or the extrapolation is not the
				// one this step would have taken. Nothing else needs restoring
				// here, because compute() always restarts from the committed
				// state
				side_commit = aux;
				// standard call
				retval = compute(true, true);
				if (retval < 0) return retval;
				implex_error = implexStressGap(stress, stress_implicit);
				// and only if the user asked for the old behaviour, fail here
				if (implex_abort_on_error && implex_error > implex_error_tolerance) {
					if (dtime_n >= implex_time_redution_limit * dtime_0) {
						retval = EC_IMPLEX_Error_Control;
					}
				}
			}
			else {
				retval = compute(true, true);
			}
		}
		else {
			retval = compute(true, true);
		}
	}

	// RECORD WHAT THIS STEP DELIVERS, in ONE place covering all four branches
	// above and not once per branch: `stress` is what getStress() is about to
	// hand the element on every path, extrapolated or implicit, and commitState()
	// is about to overwrite it with the implicit re-solve. Written per branch it
	// WAS missed, in the numeric-tangent branch - the one that is easy to forget,
	// because there only the tangent is numerical while the stress is the
	// ordinary implicit answer. On the implicit paths there is no extrapolation
	// and this equals `stress`, so the response means the same thing in both
	// regimes
	stress_implex = stress;

	// done
	return retval;
}

// ---------------------------------------------------------------------------
// get state
// ---------------------------------------------------------------------------

double ASDHysteretic1DMaterial::getStrain(void)
{
	return strain;
}

double ASDHysteretic1DMaterial::getStress(void)
{
	return stress;
}

double ASDHysteretic1DMaterial::getTangent(void)
{
	return C;
}

double ASDHysteretic1DMaterial::getInitialTangent(void)
{
	return E;
}

double ASDHysteretic1DMaterial::getEnergy(void)
{
	return energy;
}

// ---------------------------------------------------------------------------
// handle state
// ---------------------------------------------------------------------------

int ASDHysteretic1DMaterial::commitState(void)
{
	// implicit stage
	if (implex) {
		// what the extrapolated step delivered is already in stress_implex,
		// recorded by setTrialStrain
		// implicit solution. Not a second answer that gets thrown away: the
		// implicit one IS what gets committed, so this same call measures the
		// error and does the first half of the commit
		compute(false, false);
		implex_error = implexStressGap(stress_implex, stress);
		GlobalParameters::instance().setMaxError(std::max(implex_error, GlobalParameters::instance().getMaxError()));
		GlobalParameters::instance().accumulateAverageError(implex_error);
	}
	// compute energy
	energy += 0.5 * (stress_commit + stress) * (strain - strain_commit);
	// store the previously committed variables for the next move from n to n-1
	xt_commit_old = xt_commit;
	xc_commit_old = xc_commit;
	// store committed variables
	xt_commit = xt;
	xc_commit = xc;
	bt_p_commit = bt_p;
	pp_x_commit = pp_x;
	pp_y_commit = pp_y;
	bt_n_commit = bt_n;
	pn_x_commit = pn_x;
	pn_y_commit = pn_y;
	xb_p_commit = xb_p;
	xb_n_commit = xb_n;
	strain_commit = strain;
	stress_commit = stress;
	stress_eff_commit = stress_eff;
	dtime_n_commit = dtime_n;
	// done
	commit_done = true;
	return 0;
}

int ASDHysteretic1DMaterial::revertToLastCommit(void)
{
	xt = xt_commit;
	xc = xc_commit;
	bt_p = bt_p_commit;
	pp_x = pp_x_commit;
	pp_y = pp_y_commit;
	bt_n = bt_n_commit;
	pn_x = pn_x_commit;
	pn_y = pn_y_commit;
	xb_p = xb_p_commit;
	xb_n = xb_n_commit;
	strain = strain_commit;
	stress = stress_commit;
	stress_eff = stress_eff_commit;
	dtime_n = dtime_n_commit;
	return 0;
}

int ASDHysteretic1DMaterial::revertToStart(void)
{
	// state variables
	xt = 0.0;
	xt_commit = 0.0;
	xt_commit_old = 0.0;
	xc = 0.0;
	xc_commit = 0.0;
	xc_commit_old = 0.0;
	bt_p = BranchType::ReloadCross;
	pp_x = 0.0;
	pp_y = 0.0;
	bt_n = BranchType::ReloadCross;
	pn_x = 0.0;
	pn_y = 0.0;
	bt_p_commit = BranchType::ReloadCross;
	pp_x_commit = 0.0;
	pp_y_commit = 0.0;
	bt_n_commit = BranchType::ReloadCross;
	pn_x_commit = 0.0;
	pn_y_commit = 0.0;
	xb_p = 0.0;
	xb_n = 0.0;
	xb_p_commit = 0.0;
	xb_n_commit = 0.0;
	side_commit = 1;
	// implex
	dtime_n = 0.0;
	dtime_n_commit = 0.0;
	dtime_0 = 0.0;
	dtime_is_user_defined = false;
	commit_done = false;
	implex_error = 0.0;
	// strain, stress and tangent
	strain = 0.0;
	strain_commit = 0.0;
	stress = 0.0;
	stress_commit = 0.0;
	stress_implex = 0.0;
	stress_eff = 0.0;
	stress_eff_commit = 0.0;
	C = E;
	// output
	dt_bar = 0.0;
	dc_bar = 0.0;
	anchor = 0.0;
	energy = 0.0;
	return 0;
}

UniaxialMaterial* ASDHysteretic1DMaterial::getCopy(void)
{
	// the default copy constructor is safe for the member variables we use, and
	// IMPLEXObject's own copy constructor registers the copy as a NEW material
	// point (see IMPLEXManager.h)
	return new ASDHysteretic1DMaterial(*this);
}

void ASDHysteretic1DMaterial::Print(OPS_Stream& s, int flag)
{
	s << "ASDHysteretic1DMaterial - Tag: " << this->getTag() << "\n";
}

// ---------------------------------------------------------------------------
// send/recv self
// ---------------------------------------------------------------------------

int ASDHysteretic1DMaterial::sendSelf(int commitTag, Channel& theChannel)
{
	int counter;

	// aux
	int dataTag = this->getDbTag();

	// integer data
	// 12, not 10: the two limit-state COUNTS go here, because the receiving side
	// has to size its double vector before it can read the thresholds out of it
	static ID idata(12);
	counter = 0;
	idata(counter++) = dataTag;
	idata(counter++) = this->getTag();
	idata(counter++) = static_cast<int>(implex);
	idata(counter++) = static_cast<int>(implex_control);
	idata(counter++) = static_cast<int>(implex_abort_on_error);
	idata(counter++) = static_cast<int>(tangent);
	idata(counter++) = static_cast<int>(auto_regularize);
	idata(counter++) = static_cast<int>(regularization_done);
	idata(counter++) = static_cast<int>(dtime_is_user_defined);
	idata(counter++) = static_cast<int>(commit_done);
	idata(counter++) = static_cast<int>(ls_t.size());
	idata(counter++) = static_cast<int>(ls_c.size());
	if (theChannel.sendID(dataTag, commitTag, idata) < 0) {
		opserr << "ASDHysteretic1DMaterial::sendSelf() - failed to send integer data\n";
		return -1;
	}

	// double data
	int nv_state = 39;
	int nv = nv_state + 2 * Shape::NumParameters +
		ht.serializationDataSize() + hc.serializationDataSize() +
		static_cast<int>(ls_t.size() + ls_c.size());
	Vector ddata(nv);
	counter = 0;
	ddata(counter++) = E;
	ddata(counter++) = eta;
	ddata(counter++) = implex_error_tolerance;
	ddata(counter++) = implex_time_redution_limit;
	ddata(counter++) = implex_alpha;
	ddata(counter++) = lch;
	ddata(counter++) = lch_ref;
	ddata(counter++) = xt;
	ddata(counter++) = xt_commit;
	ddata(counter++) = xt_commit_old;
	ddata(counter++) = xc;
	ddata(counter++) = xc_commit;
	ddata(counter++) = xc_commit_old;
	ddata(counter++) = static_cast<double>(static_cast<int>(bt_p));
	ddata(counter++) = pp_x;
	ddata(counter++) = pp_y;
	ddata(counter++) = static_cast<double>(static_cast<int>(bt_n));
	ddata(counter++) = pn_x;
	ddata(counter++) = pn_y;
	ddata(counter++) = static_cast<double>(static_cast<int>(bt_p_commit));
	ddata(counter++) = pp_x_commit;
	ddata(counter++) = pp_y_commit;
	ddata(counter++) = static_cast<double>(static_cast<int>(bt_n_commit));
	ddata(counter++) = pn_x_commit;
	ddata(counter++) = pn_y_commit;
	ddata(counter++) = xb_p;
	ddata(counter++) = xb_n;
	ddata(counter++) = xb_p_commit;
	ddata(counter++) = xb_n_commit;
	ddata(counter++) = static_cast<double>(side_commit);
	ddata(counter++) = dtime_n;
	ddata(counter++) = dtime_n_commit;
	ddata(counter++) = dtime_0;
	ddata(counter++) = implex_error;
	ddata(counter++) = strain;
	ddata(counter++) = strain_commit;
	ddata(counter++) = stress;
	ddata(counter++) = stress_commit;
	ddata(counter++) = stress_implex;
	// nv_state ends here: stress_eff and the rest are recomputed or output only
	shape_t.serialize(ddata, counter);
	shape_c.serialize(ddata, counter);
	ht.serialize(ddata, counter);
	hc.serialize(ddata, counter);
	// the limit states LAST, so their variable count cannot shift anything else
	for (std::size_t i = 0; i < ls_t.size(); ++i)
		ddata(counter++) = ls_t[i];
	for (std::size_t i = 0; i < ls_c.size(); ++i)
		ddata(counter++) = ls_c[i];
	if (theChannel.sendVector(dataTag, commitTag, ddata) < 0) {
		opserr << "ASDHysteretic1DMaterial::sendSelf() - failed to send double data\n";
		return -1;
	}

	// done
	return 0;
}

int ASDHysteretic1DMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int counter;

	// aux
	int dataTag = this->getDbTag();

	// integer data
	static ID idata(12);
	if (theChannel.recvID(dataTag, commitTag, idata) < 0) {
		opserr << "ASDHysteretic1DMaterial::recvSelf() - failed to receive integer data\n";
		return -1;
	}
	counter = 0;
	this->setDbTag(idata(counter++));
	this->setTag(idata(counter++));
	implex = static_cast<bool>(idata(counter++));
	implex_control = static_cast<bool>(idata(counter++));
	implex_abort_on_error = static_cast<bool>(idata(counter++));
	tangent = static_cast<bool>(idata(counter++));
	auto_regularize = static_cast<bool>(idata(counter++));
	regularization_done = static_cast<bool>(idata(counter++));
	dtime_is_user_defined = static_cast<bool>(idata(counter++));
	commit_done = static_cast<bool>(idata(counter++));
	int n_ls_t = idata(counter++);
	int n_ls_c = idata(counter++);

	// double data
	int nv_state = 39;
	int nv = nv_state + 2 * Shape::NumParameters +
		ht.serializationDataSize() + hc.serializationDataSize() +
		n_ls_t + n_ls_c;
	Vector ddata(nv);
	if (theChannel.recvVector(dataTag, commitTag, ddata) < 0) {
		opserr << "ASDHysteretic1DMaterial::recvSelf() - failed to receive double data\n";
		return -1;
	}
	counter = 0;
	E = ddata(counter++);
	eta = ddata(counter++);
	implex_error_tolerance = ddata(counter++);
	implex_time_redution_limit = ddata(counter++);
	implex_alpha = ddata(counter++);
	lch = ddata(counter++);
	lch_ref = ddata(counter++);
	xt = ddata(counter++);
	xt_commit = ddata(counter++);
	xt_commit_old = ddata(counter++);
	xc = ddata(counter++);
	xc_commit = ddata(counter++);
	xc_commit_old = ddata(counter++);
	bt_p = static_cast<BranchType>(static_cast<int>(ddata(counter++)));
	pp_x = ddata(counter++);
	pp_y = ddata(counter++);
	bt_n = static_cast<BranchType>(static_cast<int>(ddata(counter++)));
	pn_x = ddata(counter++);
	pn_y = ddata(counter++);
	bt_p_commit = static_cast<BranchType>(static_cast<int>(ddata(counter++)));
	pp_x_commit = ddata(counter++);
	pp_y_commit = ddata(counter++);
	bt_n_commit = static_cast<BranchType>(static_cast<int>(ddata(counter++)));
	pn_x_commit = ddata(counter++);
	pn_y_commit = ddata(counter++);
	xb_p = ddata(counter++);
	xb_n = ddata(counter++);
	xb_p_commit = ddata(counter++);
	xb_n_commit = ddata(counter++);
	side_commit = static_cast<int>(ddata(counter++));
	dtime_n = ddata(counter++);
	dtime_n_commit = ddata(counter++);
	dtime_0 = ddata(counter++);
	implex_error = ddata(counter++);
	strain = ddata(counter++);
	strain_commit = ddata(counter++);
	stress = ddata(counter++);
	stress_commit = ddata(counter++);
	stress_implex = ddata(counter++);
	shape_t.deserialize(ddata, counter);
	shape_c.deserialize(ddata, counter);
	ht.deserialize(ddata, counter);
	hc.deserialize(ddata, counter);
	ls_t.resize(static_cast<std::size_t>(n_ls_t));
	for (int i = 0; i < n_ls_t; ++i)
		ls_t[static_cast<std::size_t>(i)] = ddata(counter++);
	ls_c.resize(static_cast<std::size_t>(n_ls_c));
	for (int i = 0; i < n_ls_c; ++i)
		ls_c[static_cast<std::size_t>(i)] = ddata(counter++);
	// with no damage the effective and the nominal stress coincide
	stress_eff = stress;
	stress_eff_commit = stress_commit;
	C = E;

	// done
	return 0;
}

// ---------------------------------------------------------------------------
// IMPL-EX error control
// ---------------------------------------------------------------------------

double ASDHysteretic1DMaterial::computeImplexErrorMetric(void)
{
	// no extrapolation, no error. This is what makes it safe for the
	// aggregation to call it on everything it has
	if (!implex)
		return 0.0;
	// the current state holds the EXPLICIT answer: the stress this step
	// delivered to the element, and the state the recorders must keep seeing
	// stress_implex already holds it - setTrialStrain records it right after the
	// extrapolation - and compute() below does not write that member, so it
	// survives the throw-away solve without being part of the peek record
	TrialState delivered;
	saveTrialState(delivered);
	// the implicit answer, at the same trial strain
	if (compute(false, false) < 0) {
		restoreTrialState(delivered);
		// no metric is not a small metric: whoever reads this must not accept
		// the step
		return std::numeric_limits<double>::quiet_NaN();
	}
	double err = implexStressGap(stress_implex, stress);   // the member
	// undo. Measuring is not allowed to move the state: the step may still be
	// rejected, and revertToLastCommit() would not put side_commit back
	restoreTrialState(delivered);
	implex_error = err;
	return err;
}

double ASDHysteretic1DMaterial::implexTimeRatio(void) const
{
	return dtime_0 > 0.0 ? dtime_n / dtime_0 : 1.0;
}

double ASDHysteretic1DMaterial::implexStressGap(double delivered, double stress_implicit) const
{
	// the reference is the largest stress this material can carry. The NORM of
	// the gap may differ from model to model (an absolute value here, the
	// largest Voigt component in 3D); the DENOMINATOR must not, or the same
	// tolerance means different things in 1D and in 3D and a stepper ends up
	// comparing numbers that are not comparable
	double ref = std::max(ht.computeMaxStress(), hc.computeMaxStress());
	if (ref <= 0.0)
		ref = 1.0;
	return std::abs(delivered - stress_implicit) / ref;
}

void ASDHysteretic1DMaterial::saveTrialState(TrialState& x) const
{
	x.xt = xt;
	x.xc = xc;
	x.bt_p = bt_p;
	x.pp_x = pp_x;
	x.pp_y = pp_y;
	x.bt_n = bt_n;
	x.pn_x = pn_x;
	x.pn_y = pn_y;
	x.xb_p = xb_p;
	x.xb_n = xb_n;
	x.side_commit = side_commit;
	x.strain = strain;
	x.stress = stress;
	x.stress_eff = stress_eff;
	x.C = C;
	x.dt_bar = dt_bar;
	x.dc_bar = dc_bar;
	x.anchor = anchor;
}

void ASDHysteretic1DMaterial::restoreTrialState(const TrialState& x)
{
	xt = x.xt;
	xc = x.xc;
	bt_p = x.bt_p;
	pp_x = x.pp_x;
	pp_y = x.pp_y;
	bt_n = x.bt_n;
	pn_x = x.pn_x;
	pn_y = x.pn_y;
	xb_p = x.xb_p;
	xb_n = x.xb_n;
	side_commit = x.side_commit;
	strain = x.strain;
	stress = x.stress;
	stress_eff = x.stress_eff;
	C = x.C;
	dt_bar = x.dt_bar;
	dc_bar = x.dc_bar;
	anchor = x.anchor;
}

// ---------------------------------------------------------------------------
// branch geometry
// ---------------------------------------------------------------------------

void ASDHysteretic1DMaterial::target(int side, double xt_in, double xc_in,
	double& tx, double& ty, double& p) const
{
	const HardeningLaw& law = (side > 0) ? ht : hc;
	double x = (side > 0) ? xt_in : xc_in;
	double x1 = (law.points().size() > 1) ? law.points()[1].x : 0.0;
	double xhat = std::max(x, x1);
	double y = law.evaluateAt(xhat).stress();
	double s = (side > 0) ? 1.0 : -1.0;
	tx = s * xhat;
	ty = s * y;
	p = s * (xhat - y / E);
}

const ASDHysteretic1DShape& ASDHysteretic1DMaterial::reloadShape(int side, BranchType bt) const
{
	// the shape of the side the excursion STARTED from: the opposite one for a
	// branch born at a zero crossing, this one for a branch born at a reversal
	int owner = (bt == BranchType::ReloadRev) ? side : -side;
	return (owner > 0) ? shape_t : shape_c;
}

void ASDHysteretic1DMaterial::unloadLeg(int side, BranchType bt, double px, double py,
	double xt_in, double xc_in, double& sx, double& sy, double& c) const
{
	double tx, ty, p_own;
	target(side, xt_in, xc_in, tx, ty, p_own);
	if (bt == BranchType::UnloadEnv) {
		sx = tx;
		sy = ty;
	}
	else {
		sx = px;
		sy = py;
	}
	// the elastic crossing of the branch's own start. It coincides with the
	// side's plastic strain unless the reversal point sits above the envelope's
	// unloading line, where p_own is the softer of the two and is taken - the
	// branch is never STIFFER than E
	double p_s = sx - sy / E;
	c = (side > 0) ? std::min(p_own, p_s) : std::max(p_own, p_s);
}

double ASDHysteretic1DMaterial::unloadSlope(int side, BranchType bt, double px, double py,
	double xt_in, double xc_in) const
{
	double sx, sy, c;
	unloadLeg(side, bt, px, py, xt_in, xc_in, sx, sy, c);
	const HardeningLaw& law = (side > 0) ? ht : hc;
	// the degenerate case returns E because that is what the polyline's
	// sentinel leg has: this reproduces eval(branch(...)) exactly, which is the
	// whole point of not building the branch
	if (std::abs(c - sx) <= law.strainTolerance())
		return E;
	return -sy / (c - sx);
}

ASDHysteretic1DMaterial::ReloadGeometry ASDHysteretic1DMaterial::reloadGeometry(
	int side, BranchType bt, double px, double py,
	double xt_in, double xc_in, bool have_k_in, double k_in) const
{
	ReloadGeometry g;
	const HardeningLaw& law = (side > 0) ? ht : hc;
	double xtol = law.strainTolerance();
	double stol = law.stressTolerance();
	double tx, ty, p_own;
	target(side, xt_in, xc_in, tx, ty, p_own);
	g.tx = tx;
	g.ty = ty;
	// the shape belongs to the side the excursion STARTED from
	int owner = (bt == BranchType::ReloadRev) ? side : -side;
	g.sh = (owner > 0) ? &shape_t : &shape_c;
	const Shape& sh = *g.sh;

	double dy_span = ty - py;
	double dy_free = dy_span;
	// FEASIBILITY - above a certain height no ABSCISSA of the vertex can
	// satisfy both hard constraints below, and they would fight each other: the
	// admissible interval shrinks with the height and closes exactly at dy_max,
	// in closed form. a4 is therefore a fraction of the height that is actually
	// AVAILABLE, not of the full stress span. Clipping dy to dy_max instead
	// would land the vertex on the single point where BOTH ceilings are tight -
	// first leg exactly k_in, second leg exactly E - which collapses the shape
	// onto the flag whatever a3 says, and makes the last reloading leg RETRACE
	// the elastic unloading line it arrives on
	if (have_k_in && k_in < E && std::abs(dy_span) > stol) {
		double room = std::abs(p_own - (px - py / E));
		double dy_max = room * k_in * E / (E - k_in);
		if (dy_max < std::abs(dy_free))
			dy_free = (dy_span > 0.0) ? dy_max : -dy_max;
	}
	double dy = sh.a4 * dy_free;
	double vy = py + dy;
	// the chord is linear in the stress, so the same fraction of the FULL span
	// (never the capped one) locates it
	double vx = px + ((std::abs(dy_span) > stol) ? (dy / dy_span) : 0.0) * (tx - px);
	double x_el_t = p_own + vy / E;   // elastic line through the TARGET
	// The peak-oriented chord can ITSELF be steeper than the incoming leg (a
	// recovered closure point moves left while the target stays put). Pull the
	// a3 = 0 REFERENCE back onto the prolongation when that happens, instead of
	// leaving it outside and letting the hard constraint below catch every a3
	// that starts from it: otherwise the whole positive half of a3 collapses
	// onto the prolongation and the response comes out with NO kink at the
	// closure point at all - the constraint is a ceiling, not a target
	if (have_k_in) {
		double x_pro = px + (vy - py) / k_in;
		vx = (side < 0) ? std::min(vx, x_pro) : std::max(vx, x_pro);
	}
	if (sh.a3 > 0.0) {
		// arrive on the target tangent to its elastic unloading line
		vx = (1.0 - sh.a3) * vx + sh.a3 * x_el_t;
	}
	else if (sh.a3 < 0.0) {
		// toward the prolongation of the incoming unloading leg: at -1 the
		// response goes through the closure point with NO change of slope. That
		// leg is the elastic one, so this is exactly "keep unloading
		// elastically past zero" - classical plasticity
		double k_ref = have_k_in ? k_in : E;
		double w = -sh.a3;
		vx = (1.0 - w) * vx + w * (px + (vy - py) / k_ref);
	}
	// HARD CONSTRAINT 1 - the kink at the closure point may never STIFFEN the
	// response: the first reloading leg is never steeper than the last
	// unloading one, otherwise the tangent jumps up while the strain keeps
	// going the same way
	if (have_k_in && std::abs(vx - px) > xtol) {
		if ((vy - py) / (vx - px) > k_in)
			vx = px + (vy - py) / k_in;
	}
	// HARD CONSTRAINT 2 - and never so flat that the SECOND leg becomes
	// stiffer than E: the vertex stays between the birth point and the elastic
	// line through the target
	double lo = (x_el_t >= px) ? px : x_el_t;
	double hi = (x_el_t >= px) ? x_el_t : px;
	if (vx < lo) vx = lo;
	if (vx > hi) vx = hi;
	if (std::abs(dy) <= stol || std::abs(vy - ty) <= stol) {
		// no vertex: the straight peak-oriented chord
		vx = px;
		vy = py;
	}
	g.vx = vx;
	g.vy = vy;
	// zero crossing: where the branch's first leg reaches zero stress
	double dx2 = vx - px;
	double dy2 = vy - py;
	if (std::abs(dx2) <= xtol || std::abs(dy2) <= stol) {
		dx2 = tx - px;
		dy2 = ty - py;
	}
	g.c = (std::abs(dx2) > xtol && std::abs(dy2) > stol) ? (px - py / (dy2 / dx2)) : px;
	return g;
}

bool ASDHysteretic1DMaterial::tipSpan(const HardeningLaw& law, const Shape& sh,
	double vx, double tx, double xtol, double& span) const
{
	if (!sh.hasTip())
		return false;
	// NO PREVIOUS PEAK, NO CUT. The virgin state is carried as a reloading
	// branch from the origin to the first backbone corner, so every reload-side
	// test lets it through: cutting that junction would round the ELASTIC LIMIT
	// on first loading and the monotonic response would stop being the
	// backbone - measured, before this guard, at 0.140 MPa on a 3.0 MPa tensile
	// strength. The target sits exactly on the first corner precisely when the
	// side has no peak beyond it to return to, which is the condition to test.
	// It also covers a reload that never got past the elastic range
	double x1 = (law.points().size() > 1) ? law.points()[1].x : 0.0;
	if (std::abs(tx) <= x1 + xtol)
		return false;
	span = tx - vx;
	return std::abs(span) > xtol;
}

int ASDHysteretic1DMaterial::tipPoints(int side, const HardeningLaw& law, const Shape& sh,
	double vx, double vy, double tx, double ty, double xtol,
	double* ox, double* oy) const
{
	double span = 0.0;
	if (!tipSpan(law, sh, vx, tx, xtol, span))
		return 0;
	// A is on the reloading leg V -> T at the fraction tip_in from the target,
	// hence also tip_in of its stress span below the backbone
	double ax = tx - sh.tip_in * span;
	double ay = ty - sh.tip_in * (ty - vy);
	// B is on the ENVELOPE ITSELF, tip_out of the same leg's strain extent past
	// the target - not on a straight prolongation of it, for two reasons: the
	// junction is continuous by construction when the engagement fires there,
	// and a backbone corner falling inside the cut cannot push the path above
	// the envelope
	double s = (side > 0) ? 1.0 : -1.0;
	double bx = tx + sh.tip_out * std::abs(span) * s;
	double by = s * law.evaluateAt(std::abs(bx)).stress();
	return roundTip(ax, ay, tx, ty, bx, by, sh.smooth, ox, oy);
}

double ASDHysteretic1DMaterial::reloadEnd(int side, BranchType bt, double px, double py,
	double xt_in, double xc_in) const
{
	// Computed in closed form, NOT by building the branch and reading its last
	// point. The cap k_in is not passed and does not need to be: it only moves
	// the vertex within the same leg, and tipSpan() is the same guard the
	// branch itself applies, so the two agree by construction rather than by
	// measurement. Building the branch here used also to need a CLAMP, because
	// a branch born on the envelope collapses into the polyline's sentinel leg
	// whose last point is c + 1: a threshold of order one, past which the
	// engagement never fires again and the response free-falls off the envelope
	// at the first peak. In closed form that case is not clamped, it is absent
	ReloadGeometry g = reloadGeometry(side, bt, px, py, xt_in, xc_in, false, 0.0);
	const HardeningLaw& law = (side > 0) ? ht : hc;
	double span = 0.0;
	double extra = 0.0;
	if (tipSpan(law, *g.sh, g.vx, g.tx, law.strainTolerance(), span))
		extra = g.sh->tip_out * std::abs(span);
	return std::abs(g.tx) + extra;
}

void ASDHysteretic1DMaterial::branch(int side, BranchType bt, double px, double py,
	double xt_in, double xc_in, bool have_k_in, double k_in, Branch& out) const
{
	const HardeningLaw& law = (side > 0) ? ht : hc;
	double xtol = law.strainTolerance();
	PointBuffer buf;

	if (!ASDHysteretic1DIsReload(bt)) {
		// unloading toward zero stress, from the envelope or from a reversal: a
		// single straight leg
		double sx, sy, c1;
		unloadLeg(side, bt, px, py, xt_in, xc_in, sx, sy, c1);
		buf.clear();
		buf.add(sx, sy);
		buf.add(c1, 0.0);
		polyline(buf, c1, xtol, E, out);
		return;
	}

	ReloadGeometry g = reloadGeometry(side, bt, px, py, xt_in, xc_in, have_k_in, k_in);
	const Shape& sh = *g.sh;
	// TIP: cut the junction with the envelope. The last leg is shortened to end
	// at A, and the path carries on past the target to B, read off the envelope,
	// so that when the engagement finally fires at B the stress is already
	// exactly the envelope's and nothing jumps
	double tipx[MaxRoundPoints], tipy[MaxRoundPoints];
	int ntip = tipPoints(side, law, sh, g.vx, g.vy, g.tx, g.ty, xtol, tipx, tipy);
	// The vertex rounding is measured against the SHORTENED leg, so the two
	// roundings meet at worst and can never overlap into each other
	double endx = (ntip > 0) ? tipx[0] : g.tx;
	double endy = (ntip > 0) ? tipy[0] : g.ty;
	double midx[MaxRoundPoints], midy[MaxRoundPoints];
	int nmid = roundVertex(px, py, g.vx, g.vy, endx, endy, sh.smooth, xtol, midx, midy);

	buf.clear();
	buf.add(px, py);
	buf.add(midx, midy, nmid);
	if (ntip > 0) {
		buf.add(tipx, tipy, ntip);
	}
	else {
		buf.add(g.tx, g.ty);
	}
	polyline(buf, g.c, xtol, E, out);
}

bool ASDHysteretic1DMaterial::neighbour(const Branch& br, double cx, double xtol,
	double& nx, double& ny) const
{
	// a branch that closes at cx, or is born there, has that point as one of
	// the two ENDS of its strain-sorted polyline; the chamfer needs the one
	// next to it, which is the direction the leg comes from or goes to
	if (std::abs(br.x[0] - cx) <= xtol) {
		nx = br.x[1];
		ny = br.y[1];
		return true;
	}
	if (std::abs(br.x[br.n - 1] - cx) <= xtol) {
		nx = br.x[br.n - 2];
		ny = br.y[br.n - 2];
		return true;
	}
	return false;
}

void ASDHysteretic1DMaterial::chamferDirect(const Shape& sh, double cx, double p_out_x,
	double xtol, double& g_in, double& g_out) const
{
	// chamfer_rs is where the cut starts, as a fraction of the stress the
	// unloading started from. The incoming leg is the ELASTIC line from that
	// point down to the corner, so a fraction of its stress is the same
	// fraction of its length: g_in is chamfer_rs with nothing to compute
	g_in = clamp01(sh.chamfer_rs);
	double span = std::abs(p_out_x - cx);
	if (span <= xtol) {
		g_out = 0.0;
		return;
	}
	// chamfer_dpe is where the cut ends, as a fraction of |DPE| past the
	// corner - and DPE IS the corner, again because the descent is elastic.
	// That is an ABSOLUTE strain offset, so unlike every other parameter it can
	// ask for more leg than there is, and is capped at the whole leg
	double want = sh.chamfer_dpe * std::abs(cx);
	double frac = want / span;
	g_out = (frac < 1.0) ? frac : 1.0;
}

void ASDHysteretic1DMaterial::chamfer(const Branch& br_in, const Branch& br_out, double cx,
	double g, double xtol, double bias, const Shape* sh,
	Branch& new_in, Branch& new_out) const
{
	double pinx, piny, poutx, pouty;
	if (!neighbour(br_in, cx, xtol, pinx, piny) ||
		!neighbour(br_out, cx, xtol, poutx, pouty)) {
		// the corner is not shared by both (the virgin state, where both sides
		// are born at the origin, is the case that matters)
		new_in = br_in;
		new_out = br_out;
		return;
	}
	double g_in, g_out;
	if (sh != 0 && (sh->chamfer_rs > 0.0 || sh->chamfer_dpe > 0.0))
		chamferDirect(*sh, cx, poutx, xtol, g_in, g_out);
	else
		Shape::chamferLegs(g, bias, g_in, g_out);
	// each cut point is a plain interpolation toward the neighbouring
	// breakpoint, so a fraction <= 1 can never overrun its own leg and no clamp
	// is needed anywhere
	double ax = cx + g_in * (pinx - cx);
	double ay = g_in * piny;
	double bx = cx + g_out * (poutx - cx);
	double by = g_out * pouty;
	if (std::abs(bx - ax) <= xtol) {
		new_in = br_in;
		new_out = br_out;
		return;
	}
	double smooth = (sh == 0) ? 0.0 : sh->smooth;
	double midx[MaxRoundPoints], midy[MaxRoundPoints];
	int nmid;
	double c;
	if (smooth > 0.0) {
		nmid = roundCorner(ax, ay, cx, bx, by, smooth, midx, midy, c);
	}
	else {
		midx[0] = ax; midy[0] = ay;
		midx[1] = bx; midy[1] = by;
		nmid = 2;
		// zero crossing of the straight cut: the new closure point
		c = (by != ay) ? (ax - ay * (bx - ax) / (by - ay)) : cx;
	}
	// everything the cut spans is replaced by it, the corner included
	double xlo = std::min(ax, bx) - xtol;
	double xhi = std::max(ax, bx) + xtol;
	PointBuffer buf;
	for (int k = 0; k < 2; ++k) {
		const Branch& src = (k == 0) ? br_in : br_out;
		Branch& dst = (k == 0) ? new_in : new_out;
		buf.clear();
		// the KEPT points first and the cut after, which is the order the
		// stable sort of polyline() resolves ties with
		for (int i = 0; i < src.n; ++i) {
			if (!(xlo <= src.x[i] && src.x[i] <= xhi))
				buf.add(src.x[i], src.y[i]);
		}
		buf.add(midx, midy, nmid);
		polyline(buf, c, xtol, E, dst);
	}
}

void ASDHysteretic1DMaterial::unloadPath(int side, BranchType bt, double px, double py,
	double xt_in, double xc_in, double& corner, Branch& cut) const
{
	Branch br;
	branch(side, bt, px, py, xt_in, xc_in, false, 0.0, br);
	corner = br.c;
	const Shape& sh = (side > 0) ? shape_t : shape_c;
	if (ASDHysteretic1DIsReload(bt) || !sh.hasChamfer()) {
		cut = br;
		return;
	}
	const HardeningLaw& law = (side > 0) ? ht : hc;
	double sigma, k;
	br.eval(br.c, br.c, sigma, k);
	Branch out;
	branch(-side, BranchType::ReloadCross, br.c, 0.0, xt_in, xc_in, true, k, out);
	Branch dummy;
	chamfer(br, out, br.c, sh.chamfer, law.strainTolerance(), sh.chamfer_bias, &sh, cut, dummy);
}

double ASDHysteretic1DMaterial::area(const Branch& br, double eps) const
{
	// exact: the branch is a polyline, so this is a sum of trapezoids over its
	// breakpoints (the ends extrapolated along the nearest leg, as everywhere)
	double lo = (br.c <= eps) ? br.c : eps;
	double hi = (br.c <= eps) ? eps : br.c;
	double total = 0.0;
	double xa = lo;
	double ya, yb, k;
	br.eval(xa, xa, ya, k);
	for (int i = 0; i < br.n; ++i) {
		double xi = br.x[i];
		if (xi <= lo || xi >= hi)
			continue;
		br.eval(xi, xi, yb, k);
		total += 0.5 * (ya + yb) * (xi - xa);
		xa = xi;
		ya = yb;
	}
	br.eval(hi, hi, yb, k);
	total += 0.5 * (ya + yb) * (hi - xa);
	return std::abs(total);
}

// ---------------------------------------------------------------------------
// the constitutive integration
// ---------------------------------------------------------------------------

int ASDHysteretic1DMaterial::compute(bool do_implex, bool do_tangent)
{
	double eps = strain;
	double x1t = (ht.points().size() > 1) ? ht.points()[1].x : 0.0;
	double x1c = (hc.points().size() > 1) ? hc.points()[1].x : 0.0;

	// get committed variables
	double xt_l = xt_commit;
	double xc_l = xc_commit;
	BranchType bt_p_l = bt_p_commit;
	double pp_x_l = pp_x_commit;
	double pp_y_l = pp_y_commit;
	BranchType bt_n_l = bt_n_commit;
	double pn_x_l = pn_x_commit;
	double pn_y_l = pn_y_commit;
	double xb_p_l = xb_p_commit;
	double xb_n_l = xb_n_commit;

	// time factor for the explicit extrapolation, CLAMPED AT ZERO. The clamp is
	// not cosmetic: the ratio is a constant-rate assumption, so a NEGATIVE one
	// extrapolates an irreversible process BACKWARDS, which is not a state this
	// material can be in. Clamping to zero degrades the step to a purely
	// elastic prediction - wrong by O(dt) like everything else here, and a
	// state that exists
	double time_factor = 1.0;
	if (implex && do_implex && (dtime_n_commit > 0.0)) {
		time_factor = dtime_n / dtime_n_commit * implex_alpha;
		if (time_factor < 0.0)
			time_factor = 0.0;
	}

	// rate coefficients (visco regularization of the measure updates)
	double rate_coeff_1 = 0.0;
	double rate_coeff_2 = 1.0;
	if ((dtime_n > 0.0) && (eta > 0.0)) {
		rate_coeff_1 = eta / (eta + dtime_n);
		rate_coeff_2 = dtime_n / (eta + dtime_n);
	}

	int side;
	double eps_leg;

	if (implex && do_implex) {
		// EXTRAPOLATED equivalent strain measures (explicit), FROZEN on/off
		// states: active side, branch types, birth points, the measure a reload
		// aims at, and the active leg (the one the committed strain sits on)
		xt_l = xt_commit + time_factor * (xt_commit - xt_commit_old);
		xc_l = xc_commit + time_factor * (xc_commit - xc_commit_old);
		side = side_commit;
		eps_leg = strain_commit;
	}
	else {
		// implicit: branch events, all decided on the COMMITTED state plus the
		// trial strain (recomputable at every trial)
		eps_leg = eps;
		double eps_n = strain_commit;
		double sig_n = stress_commit;
		double d = eps - eps_n;
		// active side from the committed stress (tie -> trial direction)
		if (sig_n > 0.0)
			side = +1;
		else if (sig_n < 0.0)
			side = -1;
		else
			side = (d >= 0.0) ? +1 : -1;

		// 1) REVERSAL event on the committed side: toward zero -> unloading
		//    branch born at the committed point (leaving it at E); away from
		//    zero -> re-aim at the side's target from it
		if (d != 0.0 && sig_n != 0.0) {
			bool away = ((sig_n > 0.0) == (d > 0.0));
			if (side > 0) {
				if (away && !ASDHysteretic1DIsReload(bt_p_l)) {
					bt_p_l = BranchType::ReloadRev;
					pp_x_l = eps_n;
					pp_y_l = sig_n;
					xb_p_l = xt_l;
				}
				else if ((!away) && ASDHysteretic1DIsReload(bt_p_l)) {
					bt_p_l = BranchType::UnloadPt;
					pp_x_l = eps_n;
					pp_y_l = sig_n;
				}
			}
			else {
				if (away && !ASDHysteretic1DIsReload(bt_n_l)) {
					bt_n_l = BranchType::ReloadRev;
					pn_x_l = eps_n;
					pn_y_l = sig_n;
					xb_n_l = xc_l;
				}
				else if ((!away) && ASDHysteretic1DIsReload(bt_n_l)) {
					bt_n_l = BranchType::UnloadPt;
					pn_x_l = eps_n;
					pn_y_l = sig_n;
				}
			}
		}

		// 1b) a committed stress of EXACTLY zero IS a crossing point: the side
		//     is the direction of travel and that side's branch is born right
		//     there. Without this the test below would measure the trial strain
		//     against the OTHER side's stale branch, find it short of ITS
		//     crossing, and re-birth the branch at that mirror point - which
		//     sits on the wrong side of the origin, so the stress jumps onto
		//     the opposite backbone. Only reachable when a step lands on
		//     sigma = 0 to the last bit, which two IDENTICAL backbones make
		//     likely: they put the two crossings at exactly opposite abscissae
		bool at_zero = (sig_n == 0.0 && d != 0.0);
		if (at_zero) {
			if (side > 0) {
				bt_p_l = BranchType::ReloadCross;
				pp_x_l = eps_n;
				pp_y_l = 0.0;
				xb_p_l = xt_l;
			}
			else {
				bt_n_l = BranchType::ReloadCross;
				pn_x_l = eps_n;
				pn_y_l = 0.0;
				xb_n_l = xc_l;
			}
		}

		// 2) ZERO CROSSING event: the other side re-aims at its target from the
		//    crossing point (the peak-oriented rule). Skipped when the
		//    committed point already IS the crossing: the branch born above
		//    starts at eps_n and the strain moves away from it, so there is no
		//    second crossing to find in this step
		if (!at_zero) {
			double corner = 0.0;
			Branch cut;
			if (side > 0) {
				unloadPath(+1, bt_p_l, pp_x_l, pp_y_l, xt_l, xc_l, corner, cut);
				if (eps < cut.c) {
					side = -1;
					bt_n_l = BranchType::ReloadCross;
					pn_x_l = corner;
					pn_y_l = 0.0;
					xb_n_l = xc_l;
				}
			}
			else {
				unloadPath(-1, bt_n_l, pn_x_l, pn_y_l, xt_l, xc_l, corner, cut);
				if (eps > cut.c) {
					side = +1;
					bt_p_l = BranchType::ReloadCross;
					pp_x_l = corner;
					pp_y_l = 0.0;
					xb_p_l = xt_l;
				}
			}
		}

		// 3) ENVELOPE ENGAGEMENT (Kuhn-Tucker in total strain) on the active
		//    side; below the first corner the measures TRACK the strain
		//    (x_pl = 0 there: the response is unaffected, but the history stays
		//    continuous for the IMPL-EX extrapolation).
		//
		//    A reloading branch with a TIP runs a little past the target before
		//    rejoining it, so the BRANCH SWITCH is held off until the end of
		//    that cut - but the MEASURE is not. The two used to be one test;
		//    the tip is the reason they had to come apart. Holding the measure
		//    back as well is the obvious reading (the path is below the
		//    envelope over the cut, so nothing seems to be damaging) and it is
		//    wrong twice over: it leaves the measure behind a strain that has
		//    already passed it, so any event that ends the reload releases the
		//    engagement at once and the stress SNAPS up onto the envelope
		//    (2.9 MPa on a 3.0 MPa material), and the catch-up is a step change
		//    in the measure, which IMPL-EX extrapolates like any other rate -
		//    the explicit error stops converging and sits at 9.5 MPa however
		//    fine the step gets. Tracking it costs nothing physically: the
		//    strain does reach past the old peak over the cut, so the damage it
		//    implies is the damage the material has actually earned
		if (side > 0) {
			if (eps > std::max(xt_l, x1t)) {
				bool in_tip =
					ASDHysteretic1DIsReload(bt_p_l) &&
					reloadShape(+1, bt_p_l).hasTip() &&
					eps <= reloadEnd(+1, bt_p_l, pp_x_l, pp_y_l, xb_p_l, xc_l);
				xt_l = rate_coeff_1 * xt_l + rate_coeff_2 * eps;
				if (!in_tip)
					bt_p_l = BranchType::UnloadEnv;
			}
			else if (eps > xt_l && xt_l < x1t) {
				// elastic zone
				xt_l = eps;
			}
		}
		else {
			if (-eps > std::max(xc_l, x1c)) {
				bool in_tip =
					ASDHysteretic1DIsReload(bt_n_l) &&
					reloadShape(-1, bt_n_l).hasTip() &&
					-eps <= reloadEnd(-1, bt_n_l, pn_x_l, pn_y_l, xt_l, xb_n_l);
				xc_l = rate_coeff_1 * xc_l + rate_coeff_2 * (-eps);
				if (!in_tip)
					bt_n_l = BranchType::UnloadEnv;
			}
			else if (-eps > xc_l && xc_l < x1c) {
				// elastic zone
				xc_l = -eps;
			}
		}
	}

	// The two branches, and the response of the active one. A branch born at a
	// zero crossing needs the slope the OTHER side's branch has there: that is
	// the reference of a negative a3 and the cap of constraint 1. The inner
	// calls are deliberately un-capped, so there is no recursion (in practice
	// the other branch is always an unloading one, which is never capped, so
	// the value is exact).
	//
	// A RELOADING branch aims at the measure it was BORN with, an unloading one
	// at the live measure. The two differ only over a tip cut, where the strain
	// has run past the target while the branch is still running into it: aiming
	// at the live measure there would make the target flee ahead of the path and
	// the cut would never close
	double own_p = ASDHysteretic1DIsReload(bt_p_l) ? xb_p_l : xt_l;
	double own_n = ASDHysteretic1DIsReload(bt_n_l) ? xb_n_l : xc_l;
	bool have_k_in_p = false, have_k_in_n = false;
	double k_in_p = 0.0, k_in_n = 0.0;
	if (bt_p_l == BranchType::ReloadCross) {
		have_k_in_p = true;
		if (ASDHysteretic1DIsReload(bt_n_l)) {
			Branch aux;
			branch(-1, bt_n_l, pn_x_l, pn_y_l, xt_l, own_n, false, 0.0, aux);
			double sigma_aux;
			aux.eval(pp_x_l, pp_x_l, sigma_aux, k_in_p);
		}
		else {
			k_in_p = unloadSlope(-1, bt_n_l, pn_x_l, pn_y_l, xt_l, own_n);
		}
	}
	if (bt_n_l == BranchType::ReloadCross) {
		have_k_in_n = true;
		if (ASDHysteretic1DIsReload(bt_p_l)) {
			Branch aux;
			branch(+1, bt_p_l, pp_x_l, pp_y_l, own_p, xc_l, false, 0.0, aux);
			double sigma_aux;
			aux.eval(pn_x_l, pn_x_l, sigma_aux, k_in_n);
		}
		else {
			k_in_n = unloadSlope(+1, bt_p_l, pp_x_l, pp_y_l, own_p, xc_l);
		}
	}
	Branch raw_p, raw_n;
	branch(+1, bt_p_l, pp_x_l, pp_y_l, own_p, xc_l, have_k_in_p, k_in_p, raw_p);
	branch(-1, bt_n_l, pn_x_l, pn_y_l, xt_l, own_n, have_k_in_n, k_in_n, raw_n);
	Branch br_p = raw_p;
	Branch br_n = raw_n;

	// CHAMFER of the closure corner, anchored to the CORNER OF THE UNLOADING
	// BRANCH - never to the opposite side's birth point, which is still the
	// previous one until the handover actually happens, and keying off it would
	// chamfer the crossing test but not the response. When the opposite side IS
	// already born there, both branches get the same A -> B segment, so the
	// response is continuous across the side change, which happens on that very
	// segment. The shape is the one of the side the excursion started from - the
	// unloading one, as for a3 and a4. Both sides can be unloading at once; the
	// two blocks are then independent and neither touches the other's branch
	for (int k = 0; k < 2; ++k) {
		int s = (k == 0) ? +1 : -1;
		const Branch& raw_own = (k == 0) ? raw_p : raw_n;
		const Branch& raw_opp = (k == 0) ? raw_n : raw_p;
		BranchType bt_own = (k == 0) ? bt_p_l : bt_n_l;
		BranchType bt_opp = (k == 0) ? bt_n_l : bt_p_l;
		double p_opp = (k == 0) ? pn_x_l : pp_x_l;
		const Shape& sh = (k == 0) ? shape_t : shape_c;
		double tol = (k == 0) ? ht.strainTolerance() : hc.strainTolerance();
		if (ASDHysteretic1DIsReload(bt_own) || !sh.hasChamfer())
			continue;
		double corner = raw_own.c;
		bool born_here = (bt_opp == BranchType::ReloadCross &&
			std::abs(p_opp - corner) <= tol);
		Branch built;
		if (!born_here) {
			double sigma_aux, k_aux;
			raw_own.eval(corner, corner, sigma_aux, k_aux);
			branch(-s, BranchType::ReloadCross, corner, 0.0, xt_l, xc_l, true, k_aux, built);
		}
		const Branch& out = born_here ? raw_opp : built;
		Branch cut_own, cut_opp;
		chamfer(raw_own, out, corner, sh.chamfer, tol, sh.chamfer_bias, &sh, cut_own, cut_opp);
		if (s > 0) {
			br_p = cut_own;
			if (born_here)
				br_n = cut_opp;
		}
		else {
			br_n = cut_own;
			if (born_here)
				br_p = cut_opp;
		}
	}

	const Branch& br = (side > 0) ? br_p : br_n;
	double sigma, k;
	br.eval(eps, eps_leg, sigma, k);

	stress = sigma;
	// with no damage the nominal and the effective stress coincide
	stress_eff = sigma;
	if (do_tangent)
		C = k;

	// stiffness reduction of the two current legs (output only)
	double sigma_aux, k_p, k_n;
	br_p.eval(eps, eps_leg, sigma_aux, k_p);
	br_n.eval(eps, eps_leg, sigma_aux, k_n);
	dt_bar = 1.0 - k_p / E;
	dc_bar = 1.0 - k_n / E;

	// store working state
	xt = xt_l;
	xc = xc_l;
	bt_p = bt_p_l;
	pp_x = pp_x_l;
	pp_y = pp_y_l;
	bt_n = bt_n_l;
	pn_x = pn_x_l;
	pn_y = pn_y_l;
	xb_p = xb_p_l;
	xb_n = xb_n_l;
	anchor = br.c;

	// save the real side if implex and !do_implex -> called from commit
	if (implex && !do_implex)
		side_commit = side;

	return 0;
}

// ---------------------------------------------------------------------------
// parameters and responses
// ---------------------------------------------------------------------------

int ASDHysteretic1DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// 1000 - elasticity & mass
	if (strcmp(argv[0], "E") == 0) {
		param.setValue(E);
		return param.addObject(1000, this);
	}
	// 2000 - lch_ref
	if (strcmp(argv[0], "lch_ref") == 0) {
		param.setValue(lch_ref);
		return param.addObject(2000, this);
	}
	// 3000 - time
	if (strcmp(argv[0], "dTime") == 0) {
		param.setValue(dtime_n);
		return param.addObject(3000, this);
	}
	if (strcmp(argv[0], "dTimeCommit") == 0) {
		param.setValue(dtime_n_commit);
		return param.addObject(3001, this);
	}
	if (strcmp(argv[0], "dTimeInitial") == 0) {
		param.setValue(dtime_0);
		return param.addObject(3002, this);
	}
	// 4000 - globals
	if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
		param.setValue(GlobalParameters::instance().getMaxError());
		return param.addObject(4000, this);
	}
	if (strcmp(argv[0], "avgImplexError") == 0 || strcmp(argv[0], "AvgImplexError") == 0) {
		param.setValue(GlobalParameters::instance().getAverageError());
		return param.addObject(4001, this);
	}
	// default
	return -1;
}

int ASDHysteretic1DMaterial::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
		// 1000 - elasticity & mass
	case 1000:
		E = info.theDouble;
		return 0;
		// 2000 - lch_ref
	case 2000:
		lch_ref = info.theDouble;
		return 0;
		// 3000 - time
	case 3000:
		dtime_n = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
	case 3001:
		dtime_n_commit = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
	case 3002:
		dtime_0 = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
		// 4000 - globals
	case 4000:
		GlobalParameters::instance().setMaxError(info.theDouble);
		return 0;
	case 4001:
		GlobalParameters::instance().setAverageError(info.theDouble);
		return 0;
	default:
		return -1;
	}
}

Vector ASDHysteretic1DMaterial::getHardeningLawVector(HardeningLawType ltype, HardeningLawPointComponent c) const
{
	Vector r;
	const HardeningLaw& law = (ltype == HardeningLawType::Tension) ? ht : hc;
	r.resize(static_cast<int>(law.points().size()));
	int counter = 0;
	for (const ASDHardeningLawPoint& p : law.points()) {
		switch (c) {
		case HardeningLawPointComponent::TotalStrain:
			r(counter++) = p.totalStrain();
			break;
		case HardeningLawPointComponent::EffectiveStress:
			r(counter++) = p.effectiveStress();
			break;
		case HardeningLawPointComponent::NominalStress:
			r(counter++) = p.stress();
			break;
		default:
			break;
		}
	}
	return r;
}

const Vector& ASDHysteretic1DMaterial::getStrainMeasure() const
{
	static Vector d(2);
	d(0) = xt;
	d(1) = xc;
	return d;
}

const Vector& ASDHysteretic1DMaterial::getEquivalentPlasticStrain() const
{
	// the ENVELOPE's plastic strain, x - y/E: with no damage, the WHOLE
	// inelastic strain. With a chamfer the path actually closes earlier, which
	// is what getClosureStrain() reports
	static Vector d(2);
	double tx, ty, p;
	target(+1, xt, xc, tx, ty, p);
	d(0) = p;
	target(-1, xt, xc, tx, ty, p);
	d(1) = -p;
	return d;
}

const Vector& ASDHysteretic1DMaterial::getClosureStrain() const
{
	// the strain at which the unloading of each side reaches zero stress: the
	// residual strain actually left behind, chamfer included. Equals the
	// equivalent plastic strain when there is no chamfer
	static Vector d(2);
	double corner;
	Branch cut;
	unloadPath(+1, BranchType::UnloadEnv, 0.0, 0.0, xt, xc, corner, cut);
	d(0) = cut.c;
	unloadPath(-1, BranchType::UnloadEnv, 0.0, 0.0, xt, xc, corner, cut);
	d(1) = -cut.c;
	return d;
}

const Vector& ASDHysteretic1DMaterial::getStiffnessReduction() const
{
	// 1 - k/E of the current leg of each side. NOT a damage variable: nothing
	// is stored and nothing is irreversible. Kept under the family's name for
	// driver and plot compatibility
	static Vector d(2);
	d(0) = dt_bar;
	d(1) = dc_bar;
	return d;
}

const Vector& ASDHysteretic1DMaterial::getImplexStress() const
{
	static Vector d(1);
	d(0) = stress_implex;
	return d;
}

namespace {

	/**
	Where a demand sits in the user's limit-state scale.

	0 = nothing reached, 1 = exactly the first limit state, 2.5 = halfway between
	the second and the third, n = the last one. A CONTINUOUS coordinate and not
	just the integer index, for two reasons: it is a smooth field to contour, and
	whoever post-processes can re-cut it at a different threshold without re-running
	the analysis - floor() of it is the discrete answer.

	WHY THE OUTPUT IS THIS SCALE AND NOT THE STRAIN. The i-th limit state sits at a
	different abscissa on each side, and at a different one again in every other
	section of the model; the index does not. So this is the quantity that can be
	compared - and maximized - across sides and across elements, which is what
	"the highest limit state ever reached" needs.

	CLAMPED AT n. Past the last threshold there is no next one to interpolate
	against, so anything above n would be an arbitrary prolongation of the last
	interval. How far past it the demand went is still readable, exactly, from the
	'equivalentTotalStrain' response.
	*/
	double ASDH1DLimitStatePosition(const std::vector<double>& ls, double x)
	{
		std::size_t n = ls.size();
		if (n == 0 || !(x > 0.0))
			return 0.0;
		if (x >= ls[n - 1])
			return static_cast<double>(n);
		// below the first one: the fraction of the way to it. Not zero, so that the
		// field is continuous and says how close the point is to the first state
		if (x < ls[0])
			return x / ls[0];
		for (std::size_t k = 1; k < n; ++k) {
			if (x < ls[k])
				return static_cast<double>(k) +
					(x - ls[k - 1]) / (ls[k] - ls[k - 1]);
		}
		return static_cast<double>(n);
	}

}

const Vector& ASDHysteretic1DMaterial::getLimitStateRatio() const
{
	// THE MAXIMUM EVER REACHED, and it needs no state variable of its own: xt and
	// xc are the equivalent strain measures, which are monotone non-decreasing and
	// are exactly the peak |strain| reached on their side. Measured on a history
	// with shrinking amplitude: largest decrease 0.000e+00 over 140 steps, and
	// xt/ey == max(+strain)/ey to the digit.
	//
	// Deriving it instead of latching a flag is also what makes it CORRECT under
	// rollback: a latched flag would keep a state the analysis rolled back with
	// revertToLastCommit, and under IMPL-EX it would keep an EXTRAPOLATED peak that
	// the implicit commit does not confirm. This follows both for free.
	//
	// Read from xt/xc and not from the committed pair, so that this response and
	// 'equivalentTotalStrain' always tell the same story.
	static Vector d(1);
	d(0) = std::max(ASDH1DLimitStatePosition(ls_t, xt),
		ASDH1DLimitStatePosition(ls_c, xc));
	return d;
}

const Vector& ASDHysteretic1DMaterial::getImplexError() const
{
	static Vector d(1);
	d(0) = implex_error;
	return d;
}

const Vector& ASDHysteretic1DMaterial::getTimeIncrements() const
{
	static Vector d(3);
	d(0) = dtime_n;
	d(1) = dtime_n_commit;
	d(2) = dtime_0;
	return d;
}

const Vector& ASDHysteretic1DMaterial::getBranchState() const
{
	static Vector d(8);
	d(0) = static_cast<double>(static_cast<int>(bt_p));
	d(1) = pp_x;
	d(2) = pp_y;
	d(3) = static_cast<double>(static_cast<int>(bt_n));
	d(4) = pn_x;
	d(5) = pn_y;
	d(6) = static_cast<double>(side_commit);
	d(7) = anchor;
	return d;
}

const Vector& ASDHysteretic1DMaterial::getFreeEnergy() const
{
	// Helmholtz free energy: the area under the path the material would ACTUALLY
	// follow from the current state down to zero stress - the remainder of the
	// current branch when that branch is already unloading, the unloading branch
	// that would be BORN here when it is reloading, since reversing is what one
	// does to unload from a reloading branch. Both are polylines, so the
	// integral is exact.
	//
	// With every shape parameter at zero each unloading branch is the elastic
	// line and this collapses to the familiar sigma^2/(2E). A non-zero chamfer
	// moves the zero-stress crossing past the plastic strain, so the descent
	// releases MORE than that and the closed form stops being the free energy:
	// what is recoverable is what the material would actually give back
	static Vector d(1);
	d(0) = 0.0;
	if (stress != 0.0) {
		int side = (stress > 0.0) ? +1 : -1;
		BranchType bt = (side > 0) ? bt_p : bt_n;
		double px = (side > 0) ? pp_x : pn_x;
		double py = (side > 0) ? pp_y : pn_y;
		if (ASDHysteretic1DIsReload(bt)) {
			bt = BranchType::UnloadPt;
			px = strain;
			py = stress;
		}
		double corner;
		Branch cut;
		unloadPath(side, bt, px, py, xt, xc, corner, cut);
		d(0) = area(cut, strain);
	}
	return d;
}

const Vector& ASDHysteretic1DMaterial::getShapeVector(int side) const
{
	static Vector d(Shape::NumParameters);
	const Shape& sh = (side > 0) ? shape_t : shape_c;
	int pos = 0;
	sh.serialize(d, pos);
	return d;
}

Response* ASDHysteretic1DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	// utils
	auto make_resp = [&output, this](int rid, const Vector& v, const std::vector<std::string>* labels = nullptr) -> MaterialResponse* {
		output.tag("UniaxialMaterialOutput");
		output.attr("matType", this->getClassType());
		output.attr("matTag", this->getTag());
		if (labels) {
			for (const auto& item : (*labels))
				output.tag("ResponseType", item.c_str());
		}
		MaterialResponse* resp = new MaterialResponse(this, rid, v);
		output.endTag();
		return resp;
	};

	// labels
	static std::vector<std::string> lb_eqpl_strain = { "PLE+", "PLE-" };
	static std::vector<std::string> lb_tot_strain = { "TE+", "TE-" };
	static std::vector<std::string> lb_closure = { "CL+", "CL-" };
	static std::vector<std::string> lb_reduction = { "d+", "d-" };
	static std::vector<std::string> lb_implex_error = { "Error" };
	static std::vector<std::string> lb_implex_stress = { "sigma11" };
	static std::vector<std::string> lb_limit_state = { "LS" };
	static std::vector<std::string> lb_time = { "dTime", "dTimeCommit", "dTimeInitial" };
	static std::vector<std::string> lb_shape = { "a3", "a4", "chamfer", "chamferBias",
		"chamferRs", "chamferDpe", "smooth", "tipIn", "tipOut" };
	static std::vector<std::string> lb_branch = { "btP", "pxP", "pyP", "btN", "pxN", "pyN", "side", "anchor" };
	static std::vector<std::string> lb_energy = { "psi" };

	// check specific responses
	if (argc > 0) {
		// 1000 - compressive backbone
		if (strcmp(argv[0], "Ce") == 0)
			return make_resp(1000, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::TotalStrain));
		if (strcmp(argv[0], "Cs") == 0)
			return make_resp(1001, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::NominalStress));
		// 1100 - tensile backbone
		if (strcmp(argv[0], "Te") == 0)
			return make_resp(1100, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::TotalStrain));
		if (strcmp(argv[0], "Ts") == 0)
			return make_resp(1101, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::NominalStress));
		// 2000 - state
		if (strcmp(argv[0], "equivalentPlasticStrain") == 0 || strcmp(argv[0], "EquivalentPlasticStrain") == 0)
			return make_resp(2001, getEquivalentPlasticStrain(), &lb_eqpl_strain);
		if (strcmp(argv[0], "equivalentTotalStrain") == 0 || strcmp(argv[0], "EquivalentTotalStrain") == 0)
			return make_resp(2002, getStrainMeasure(), &lb_tot_strain);
		if (strcmp(argv[0], "closureStrain") == 0 || strcmp(argv[0], "ClosureStrain") == 0)
			return make_resp(2006, getClosureStrain(), &lb_closure);
		// 'damage' is the family's name for what here is only the instantaneous
		// stiffness reduction of the active leg: no state is attached to it
		if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0 ||
			strcmp(argv[0], "stiffnessReduction") == 0)
			return make_resp(2000, getStiffnessReduction(), &lb_reduction);
		// 2100 - the branch state machine
		if (strcmp(argv[0], "branch") == 0 || strcmp(argv[0], "Branch") == 0)
			return make_resp(2100, getBranchState(), &lb_branch);
		// 2200 - the shape of each side
		if (strcmp(argv[0], "shapeT") == 0 || strcmp(argv[0], "ShapeT") == 0)
			return make_resp(2200, getShapeVector(+1), &lb_shape);
		if (strcmp(argv[0], "shapeC") == 0 || strcmp(argv[0], "ShapeC") == 0)
			return make_resp(2201, getShapeVector(-1), &lb_shape);
		// 2300 - the recoverable energy
		if (strcmp(argv[0], "freeEnergy") == 0 || strcmp(argv[0], "FreeEnergy") == 0)
			return make_resp(2300, getFreeEnergy(), &lb_energy);
		// 3000 - implex error
		if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0)
			return make_resp(3000, getImplexError(), &lb_implex_error);
		// the stress the step DELIVERED, which after the commit is otherwise
		// gone: under IMPL-EX the implicit re-solve installs its own over it
		if (strcmp(argv[0], "implexStress") == 0 || strcmp(argv[0], "ImplexStress") == 0)
			return make_resp(3003, getImplexStress(), &lb_implex_stress);
		// 5000 - THE HIGHEST USER LIMIT STATE EVER REACHED, as a continuous
		// position in the limit-state scale: see getLimitStateRatio()
		if (strcmp(argv[0], "limitStateRatio") == 0 || strcmp(argv[0], "LimitStateRatio") == 0 ||
			strcmp(argv[0], "LS") == 0) {
			if (ls_t.size() == 0 && ls_c.size() == 0) {
				// SAY IT, rather than answer a zero forever. A recorder that gets a
				// valid response full of zeros looks like a model that never yielded
				opserr << "uniaxialMaterial ASDHysteretic1D Warning: the "
					"'limitStateRatio' response was asked of material " << getTag()
					<< ", which was given no limit states. Pass them with "
					"'-limitStates $x1 $x2 ...', ascending, in STRAIN units.\n";
				return nullptr;
			}
			return make_resp(5000, getLimitStateRatio(), &lb_limit_state);
		}
		// 4000 - internal time
		if (strcmp(argv[0], "time") == 0 || strcmp(argv[0], "Time") == 0)
			return make_resp(4000, getTimeIncrements(), &lb_time);
	}

	// otherwise return base-class response
	return UniaxialMaterial::setResponse(argv, argc, output);
}

int ASDHysteretic1DMaterial::getResponse(int responseID, Information& matInformation)
{
	switch (responseID) {
	case 1000: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::TotalStrain));
	case 1001: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::NominalStress));
	case 1100: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::TotalStrain));
	case 1101: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::NominalStress));
	case 2000: return matInformation.setVector(getStiffnessReduction());
	case 2001: return matInformation.setVector(getEquivalentPlasticStrain());
	case 2002: return matInformation.setVector(getStrainMeasure());
	case 2006: return matInformation.setVector(getClosureStrain());
	case 2200: return matInformation.setVector(getShapeVector(+1));
	case 2201: return matInformation.setVector(getShapeVector(-1));
	case 3000: return matInformation.setVector(getImplexError());
	case 3003: return matInformation.setVector(getImplexStress());
	case 5000: return matInformation.setVector(getLimitStateRatio());
	case 4000: return matInformation.setVector(getTimeIncrements());
	case 2100: return matInformation.setVector(getBranchState());
	case 2300: return matInformation.setVector(getFreeEnergy());
	default:
		break;
	}
	// THE BASE CLASS OWNS stress/strain/tangent/energy, and setResponse above
	// already hands them out by falling through to it: refusing them here
	// instead of forwarding leaves a recorder holding a valid Response that is
	// never filled, so it records nothing and says nothing
	return UniaxialMaterial::getResponse(responseID, matInformation);
}
