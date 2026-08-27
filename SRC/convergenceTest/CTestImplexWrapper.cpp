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

#include <CTestImplexWrapper.h>
#include <IMPLEXManager.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <EquiSolnAlgo.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <ParallelAgreement.h>

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#include <mpi.h>
#endif
#if defined(_PARALLEL_PROCESSING)
// the PartitionedDomain flag, and ONLY the build that has a PartitionedDomain
// may read it - see aggregateImplexError()
extern bool OPS_PARTITIONED;
#endif

CTestImplexWrapper::CTestImplexWrapper()
	: ConvergenceTest(CONVERGENCE_TEST_CTestImplexWrapper)
	, norms(1)
{
}

CTestImplexWrapper::CTestImplexWrapper(ConvergenceTest* test, double maxError,
	double maxFrac, double maxReduction, FloorPolicy floor, int print)
	: ConvergenceTest(CONVERGENCE_TEST_CTestImplexWrapper)
	, theTest(test)
	, maxImplexError(maxError)
	, maxFraction(maxFrac)
	, maxReductionFactor(maxReduction)
	, onFloor(floor)
	, printFlag(print)
	, norms(1)
{
}

CTestImplexWrapper::~CTestImplexWrapper()
{
	if (theTest)
		delete theTest;
}

ConvergenceTest* CTestImplexWrapper::getCopy(int iterations)
{
	// the inner test is copied too, not shared: two wrappers that share one
	// test would reset each other's iteration count on start()
	// EVERY option travels, and maxFraction is the one that must: Broyden, BFGS
	// and NewtonLineSearch keep a second test from getCopy(), so a copy that
	// dropped it would run a different criterion - and, if the reduction below
	// were ever gated on it, a different number of MPI collectives inside one
	// process
	ConvergenceTest* innerCopy = theTest ? theTest->getCopy(iterations) : 0;
	return new CTestImplexWrapper(innerCopy, maxImplexError, maxFraction,
		maxReductionFactor, onFloor, printFlag);
}

int CTestImplexWrapper::setEquiSolnAlgo(EquiSolnAlgo& theAlgo)
{
	if (theTest == 0) {
		opserr << "WARNING: CTestImplexWrapper::setEquiSolnAlgo() - no inner test set.\n";
		return -1;
	}
	return theTest->setEquiSolnAlgo(theAlgo);
}

int CTestImplexWrapper::start(void)
{
	if (theTest == 0) {
		opserr << "WARNING: CTestImplexWrapper::start() - no inner test set.\n";
		return -1;
	}
	// a new attempt at a step: forget which material points took part in the
	// previous one. This is called once per solveCurrentStep by every
	// equiSolnAlgo, which is exactly the scope we want
	IMPLEXManager::instance().clearTouched();
	lastImplexError = 0.0;
	return theTest->start();
}

CTestImplexWrapper::Measure CTestImplexWrapper::aggregateImplexError(void)
{
	Measure m;

	// measure this rank's material points. Only the ones that took part in
	// this step, and only those not measured yet: a second call in the same
	// step costs nothing.
	//
	// The tolerance goes IN so that the count of violations comes out of the
	// one pass that is being paid for anyway. The registry does not decide
	// anything with it - it compares.
	const IMPLEXManager::Aggregate& agg =
		IMPLEXManager::instance().aggregate(maxImplexError);
	double e = agg.any_nan ? std::numeric_limits<double>::quiet_NaN() : agg.max;
	// as doubles from here on: they are about to cross MPI, where one datatype
	// is better than two, and a count is exact in a double up to 2^53
	m.over = static_cast<double>(agg.count_over);
	m.count = static_cast<double>(agg.count);

	// HOW FAR THE SCHEDULE HAS CUT, measured against THIS analysis step's
	// nominal - see nominalDt in the header for why the materials' own ratio
	// cannot answer this once a second analysis step exists.
	//
	// ops_Dt is the same quantity the materials extrapolate with (they take
	// dtime_n = ops_Dt), and Domain::update() has already set it for the
	// current step by the time any convergence test runs. When it is not
	// usable - no integrator has set a pseudo-time increment, or a driver has
	// pinned the material's dt through a Parameter - fall back to what the
	// registry reports, which is the previous behaviour.
	double dtNow = ops_Dt;
	if (dtNow > nominalDt)
		nominalDt = dtNow;
	double minTimeRatio = (dtNow > 0.0 && nominalDt > 0.0) ?
		(dtNow / nominalDt) : agg.min_time_ratio;

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
	// THE REDUCTION LIVES HERE, not in the registry, because this is the only
	// place that knows the call is collective: test() is reached by every rank
	// the same number of times, while the registry is just an accumulator with
	// no business knowing about MPI.
	//
	// WHO IS EXEMPT, and it is NOT the same in the two parallel builds - reading
	// OPS_PARTITIONED in the wrong one is what used to leave rank 0 out of an
	// Allreduce that every other rank had entered. MPI does not report that as a
	// missing participant: it pairs rank 0's NEXT collective on the communicator
	// - worstStepResult()'s MPI_2INT/MPI_MINLOC - with this MPI_DOUBLE/MPI_MAX
	// and aborts on the datatype mismatch, naming the innocent call.
	//
	//  * OpenSeesMP (_PARALLEL_INTERPRETERS): every rank runs the script and
	//    builds its own copy of this test, so every rank gets here and the
	//    reduction is UNCONDITIONAL. OPS_PARTITIONED is a PartitionedDomain flag
	//    that nothing in this build ever sets - partitionModel() is compiled out
	//    - so it is not just false, it is meaningless here.
	//  * OpenSeesSP (_PARALLEL_PROCESSING): before the model is partitioned the
	//    object exists on P0 only and a collective would hang, so P0 alone stays
	//    out until then. Once partitioned every rank has its own copy -
	//    specifyCTest sends it to the subdomains - so all of them get here.
	//
	// This is the shape AutoConstraintHandler uses for the same situation.
	//
	// PRECONDITION, and it belongs to the caller: test() must be reached the same
	// number of times on every rank. It holds when the inner test's verdict is
	// global, which is what the stock tests measure - MumpsParallelSOE
	// broadcasts X and Allreduces B, so every rank sees the same norms - and the
	// parallel Newton already depends on it. Three ways to break it, none of them
	// detectable from in here: an inner test that measures rank-local data (one
	// written in-house, CTestPFEM), the "N independent analyses, one per rank"
	// use of OpenSeesMP, and the legacy -implexAbort, where a material returning
	// EC_IMPLEX_Error_Control on ONE rank fails its element's update and that
	// rank leaves the algorithm before ever reaching test().
	//
	// NOTE for OpenSeesSP: this also assumes master and subdomain analyses reach
	// test() in lockstep. If they do not, this hangs - loudly, and never with
	// a silently wrong answer, which is why it is worth trying rather than
	// refusing to run. A partial maximum would be the bad outcome: ranks
	// rejecting different steps.
	int pid = 0, np = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	bool do_allreduce = true;
	if (np == 1) do_allreduce = false;
#if defined(_PARALLEL_PROCESSING)
	if (pid == 0 && !OPS_PARTITIONED) do_allreduce = false;
#endif // defined(_PARALLEL_PROCESSING)
#if defined(_PARALLEL_INTERPRETERS)
	// and nobody is exempt in OpenSeesMP, but a run whose processes do not share a
	// model owes no collective at all: that is the "N independent analyses, one per
	// rank" case named above, where reducing would pair this test's maximum with an
	// unrelated analysis's - see ParallelAgreement.h
	if (!OPS_inCoupledParallelRun()) do_allreduce = false;
#endif // defined(_PARALLEL_INTERPRETERS)
	if (do_allreduce) {
		double local[2] = { std::isnan(e) ? 1.0e300 : e, -minTimeRatio };
		double global[2] = { 0.0, 0.0 };
		if (MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) == MPI_SUCCESS) {
			// 1.0e300 travels in place of a NaN, which MPI_MAX does not order
			e = (global[0] >= 1.0e300) ? std::numeric_limits<double>::quiet_NaN() : global[0];
			minTimeRatio = -global[1];
		}
		else {
			opserr << "CTestImplexWrapper: MPI_Allreduce failed, the IMPL-EX error is this rank's only\n";
		}

		// THE COUNTS NEED A SUM, AND A SUM IS NOT A MAX. So they take a second
		// collective, and three things about it are not free choices:
		//
		//  * it is INSIDE the same guard, right after the first, so the number
		//    of collectives a rank performs depends on nothing but that guard.
		//    A collective whose existence depended on a policy value - even one
		//    as rank-uniform as maxFraction - is the exact shape of the bug
		//    that left rank 0 out of the reduction above;
		//  * there is no early return between the two, not even when the first
		//    one fails. AutoConstraintHandler does four in a row this way, warn
		//    and carry on, and that is the model. A return here would be a rank
		//    leaving with a collective still owed;
		//  * it runs even when e is NaN. isnan(e) is a RANK-LOCAL fact before
		//    the reduction, so deciding anything on it before both collectives
		//    are done is how the ranks stop agreeing.
		//
		// It costs two doubles per measured step - once per attempted step, not
		// per iteration - and it changes no number when maxFraction is 0.
		double lcount[2] = { m.over, m.count };
		double gcount[2] = { 0.0, 0.0 };
		if (MPI_Allreduce(lcount, gcount, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) == MPI_SUCCESS) {
			m.over = gcount[0];
			m.count = gcount[1];
		}
		else {
			opserr << "CTestImplexWrapper: MPI_Allreduce failed, the IMPL-EX point count is this rank's only\n";
		}
	}
#endif

	// after the reduction, never before: a fraction of a rank's own points is
	// not the fraction of the model's
	m.max = e;
	m.minTimeRatio = minTimeRatio;
	m.fraction = (m.count > 0.0) ? (m.over / m.count) : 0.0;
	return m;
}

int CTestImplexWrapper::test(void)
{
	if (theTest == 0) {
		opserr << "WARNING: CTestImplexWrapper::test() - no inner test set.\n";
		return -2;
	}

	// the test the user chose comes first, and it can already fail or ask for
	// another iteration
	int result = theTest->test();
	if (result == -1)
		return -1;                       // not equilibrated: no metric to take
	if (result < 0)
		return result;                   // it failed on its own terms

	// converged: now, and only now, the IMPL-EX criterion
	Measure m = aggregateImplexError();
	double e = m.max;
	lastImplexError = e;

	if (std::isnan(e)) {
		// no metric is not a small metric: a material that could not solve
		// implicitly cannot say anything about its own error
		if (printFlag != 0)
			opserr << "CTestImplexWrapper: no IMPL-EX metric (NaN), step rejected\n";
		return -2;
	}

	// HOW MANY points are over tolerance, not how bad the worst one is. At
	// maxFraction 0 the two are the same question - the fraction is nonzero
	// exactly when the maximum is over tolerance, and 'exactly' is not an
	// approximation: the maximum IS one of the points that were counted, so
	// there is no comparison here that the old 'e <= maxImplexError' answered
	// differently. An empty population gives a fraction of zero and is
	// accepted, as it always was
	if (m.fraction <= maxFraction)
		return result;                   // both criteria satisfied

	if (onFloor == Floor_Accept && m.minTimeRatio <= maxReductionFactor) {
		// the step cannot be reduced any further: taking it is better than
		// deadlocking. The other policy, which is the one to use under an
		// external adaptive scheme, is to fail and let it decide
		if (printFlag != 0)
			opserr << "CTestImplexWrapper: IMPL-EX error " << e << " > " << maxImplexError
			<< " at " << m.over << " of " << m.count << " points"
			<< " but the step is at the floor (" << m.minTimeRatio << "), accepted\n";
		return result;
	}

	if (printFlag != 0) {
		// the counts, not just the worst value: 'over at 1 of 240000 points' and
		// 'over at 900 of 240000' are the same fraction to nobody, and they are
		// the two cases maxFraction exists to tell apart
		opserr << "CTestImplexWrapper: IMPL-EX error " << e << " > " << maxImplexError
			<< " at " << m.over << " of " << m.count << " points";
		if (maxFraction > 0.0)
			opserr << " (fraction " << m.fraction << " > " << maxFraction << ")";
		opserr << ", step rejected\n";
	}

	// -2 and never -1: iterating does not fix a discretisation error
	return -2;
}

int CTestImplexWrapper::getNumTests(void)
{
	return theTest ? theTest->getNumTests() : 0;
}

int CTestImplexWrapper::getMaxNumTests(void)
{
	return theTest ? theTest->getMaxNumTests() : 0;
}

double CTestImplexWrapper::getRatioNumToMax(void)
{
	return theTest ? theTest->getRatioNumToMax() : 0.0;
}

const Vector& CTestImplexWrapper::getNorms(void)
{
	if (theTest == 0) {
		norms.resize(1);
		norms(0) = lastImplexError;
		return norms;
	}
	const Vector& inner = theTest->getNorms();
	int n = inner.Size();
	if (norms.Size() != n + 1)
		norms.resize(n + 1);
	for (int i = 0; i < n; ++i)
		norms(i) = inner(i);
	norms(n) = lastImplexError;
	return norms;
}

int CTestImplexWrapper::sendSelf(int cTag, Channel& theChannel)
{
	// own data, plus the class tag of the inner test so that the broker on the
	// other side can build one to receive into.
	//
	// maxFraction travels because it is POLICY, and a subdomain running a
	// different criterion from its master would reject different steps - the
	// one failure mode this whole reduction exists to prevent. nominalDt does
	// not travel, because it is state: see the header
	static Vector x(6);
	x(0) = maxImplexError;
	x(1) = maxFraction;
	x(2) = maxReductionFactor;
	x(3) = static_cast<double>(static_cast<int>(onFloor));
	x(4) = static_cast<double>(printFlag);
	x(5) = static_cast<double>(theTest ? theTest->getClassTag() : -1);
	if (theChannel.sendVector(this->getDbTag(), cTag, x) < 0) {
		opserr << "CTestImplexWrapper::sendSelf() - failed to send data\n";
		return -1;
	}
	if (theTest) {
		if (theTest->sendSelf(cTag, theChannel) < 0) {
			opserr << "CTestImplexWrapper::sendSelf() - failed to send the inner test\n";
			return -1;
		}
	}
	return 0;
}

int CTestImplexWrapper::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	static Vector x(6);
	if (theChannel.recvVector(this->getDbTag(), cTag, x) < 0) {
		opserr << "CTestImplexWrapper::recvSelf() - failed to receive data\n";
		return -1;
	}
	maxImplexError = x(0);
	maxFraction = x(1);
	maxReductionFactor = x(2);
	onFloor = static_cast<FloorPolicy>(static_cast<int>(x(3)));
	printFlag = static_cast<int>(x(4));
	int innerClassTag = static_cast<int>(x(5));
	// rebuild the inner test only if what we have is not already of the right
	// type, as the domain-decomposition analyses do with theirs
	if (theTest != 0 && theTest->getClassTag() != innerClassTag) {
		delete theTest;
		theTest = 0;
	}
	if (theTest == 0 && innerClassTag >= 0) {
		theTest = theBroker.getNewConvergenceTest(innerClassTag);
		if (theTest == 0) {
			opserr << "CTestImplexWrapper::recvSelf() - the broker could not build a test of class "
				<< innerClassTag << "\n";
			return -1;
		}
	}
	if (theTest != 0) {
		if (theTest->recvSelf(cTag, theChannel, theBroker) < 0) {
			opserr << "CTestImplexWrapper::recvSelf() - failed to receive the inner test\n";
			return -1;
		}
	}
	return 0;
}
