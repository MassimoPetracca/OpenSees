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
// The registry of the objects that integrate with IMPL-EX, and the one place
// where their error metrics are aggregated.
//
// WHY A REGISTRY AND NOT A WALK OVER THE DOMAIN. Whoever wants the IMPL-EX
// error of a model - today the convergence test wrapper - has to reach every
// material point. Going through setResponse means guessing an addressing that
// is element-specific and does not even agree on the index base ('section i'
// is 1-based on the beams, 'fiber j' is 0-based on the fiber sections), with
// no way to know how many points there are: an element that does not forward
// the key DISAPPEARS from the aggregate in silence, and an aggregate that
// loses the wrong points is worse than no aggregate. Objects that carry an
// IMPL-EX scheme know they do; they say so here, once, in their constructor.
//
// THREADING. OpenSees has no shared-memory parallelism inside a rank, so this
// needs no locking. Across ranks the aggregate is partial by construction and
// the reduction belongs to the consumer, which is the only one that knows the
// call is collective at that point - see CTestImplexWrapper.
//
// WHAT 'touched' IS FOR. Only the objects that were actually updated in the
// current step take part in the aggregate. That excludes, with one mechanism
// and no special cases: the prototypes held by the model builder, which are
// never assigned to a gauss point; the objects of DEACTIVATED elements, whose
// update is skipped by staged construction and whose frozen error would
// otherwise keep the maximum pinned to a stale value forever; and anything
// else that skips its update for reasons of its own.

#ifndef IMPLEXManager_h
#define IMPLEXManager_h

#include <vector>
#include <cstddef>

class IMPLEXManager;

/**
The interface an object exposes to take part in the IMPL-EX error control.
Materials and elements alike: nothing here says 'material'.
*/
class IMPLEXObject {
	friend class IMPLEXManager;

public:
	IMPLEXObject();
	// a copy is a NEW object with its own slot: identity is never copied, or
	// two objects would share one entry of the registry and the second one to
	// die would remove somebody else's
	IMPLEXObject(const IMPLEXObject&);
	IMPLEXObject& operator = (const IMPLEXObject&) { return *this; }
	virtual ~IMPLEXObject();

public:
	/**
	Measure the extrapolation error at the CURRENT trial state and return it,
	normalized, WITHOUT changing the state.

	Non-destructive is not a nicety. It is called between convergence and
	commit, so (a) if the step is rejected the object must be exactly as it
	was - and revertToLastCommit() cannot help, because it does not restore
	the quantities that are FROZEN AT COMMIT and that an implicit pass
	overwrites - and (b) the recorders run right after, and must see the state
	the step delivered.

	Returns 0 for an object that is not running IMPL-EX, which is what makes
	it safe to call unconditionally.
	*/
	virtual double computeImplexErrorMetric() = 0;

	/**
	The size of the current step over the size of the first one, which is what
	a policy compares against its floor. 1 when the object does not know.
	*/
	virtual double implexTimeRatio() const { return 1.0; }

protected:
	// say that this object took part in the current step. One store, called
	// from setTrialStrain (a material) or update (an element)
	inline void implexTouch() { m_implex_touched = true; }

private:
	std::size_t m_implex_uid = 0;
	bool m_implex_touched = false;
};

/**
The registry. One instance per process, which is the right scope: each rank
owns its own objects.
*/
class IMPLEXManager {
public:
	/**
	What an aggregation pass found. The consumer owns the policy: which of
	these to use, against which tolerance, and what to do at the floor.
	*/
	struct Aggregate {
		// the largest error over the objects that took part
		double max = 0.0;
		// HOW MANY of them are strictly over the threshold the CONSUMER handed
		// to aggregate(). The registry does not interpret that number and owns
		// no policy: it counts.
		//
		// This is the other half of the criterion, and a COUNT rather than a
		// mean on purpose. The metric is bimodal - exactly zero on the elastic
		// steps, O(tol) on the few points that carry plastic flow - so a mean
		// over material points is diluted by the elastic ones AND depends on
		// the mesh: with N points and k active it is (k/N) times the mean over
		// the active ones, so the same physics read on a finer mesh gives a
		// smaller number and a tolerance means something different on every
		// model. count_over/count is a fraction: dimensionless, independent of
		// the model's size, and zero exactly when the maximum is within
		// tolerance - which is how the plain 'reject if the worst point is
		// over' criterion is the fraction criterion at threshold zero.
		std::size_t count_over = 0;
		// sum and count. count is also the denominator of the fraction above;
		// sum is for an average as a DIAGNOSTIC only. Note that the iteration
		// order of the registry is not stable, so a sum is not bit-reproducible;
		// a maximum and a count are
		double sum = 0.0;
		std::size_t count = 0;
		// a metric that is not a number is never 'within tolerance'
		bool any_nan = false;
		// the smallest step ratio among the objects that took part
		double min_time_ratio = 1.0;
	};

private:
	IMPLEXManager() = default;
	IMPLEXManager(const IMPLEXManager&) = delete;
	IMPLEXManager& operator = (const IMPLEXManager&) = delete;

public:
	static IMPLEXManager& instance();

	/**
	Start of a step: forget who took part and reset the accumulated numbers.
	*/
	void clearTouched();

	/**
	Measure the objects that took part SINCE THE LAST CALL, accumulate, and
	forget them.

	`overThreshold` is counted against, not interpreted: count_over ends up
	holding how many of the measured objects came out strictly above it. Pass
	the tolerance to get the fraction the convergence test needs; pass zero to
	get how many are active at all.

	Measuring only the newly touched ones is what makes a second call in the
	same step free: an algorithm that asks twice without any new
	setTrialStrain in between gets the same numbers back without a single
	extra constitutive pass. Broyden, BFGS and NewtonLineSearch are the ones
	that ask twice, and not from two call sites of their own - they build a
	SECOND convergence test with getCopy(), which for a CTestImplexWrapper is
	a second wrapper sharing this one registry. In the Broyden/BFGS chain the
	second ask is free because nothing updates the domain between the two.
	*/
	const Aggregate& aggregate(double overThreshold);

	/**
	The last aggregate, without measuring anything.
	*/
	inline const Aggregate& lastAggregate() const { return m_aggregate; }

	inline std::size_t size() const { return m_objects.size(); }

private:
	std::size_t add(IMPLEXObject* obj);
	void remove(IMPLEXObject* obj);
	friend class IMPLEXObject;

private:
	// slot vector, not a map: one pointer per object instead of a node per
	// object, and a contiguous walk. Removal is a swap with the last one, so
	// the moved object is told its new slot
	std::vector<IMPLEXObject*> m_objects;
	Aggregate m_aggregate;
};

#endif // IMPLEXManager_h
