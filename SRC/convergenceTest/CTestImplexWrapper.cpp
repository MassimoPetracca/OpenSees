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

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
extern bool OPS_PARTITIONED;
#include <mpi.h>
#endif

CTestImplexWrapper::CTestImplexWrapper()
	: ConvergenceTest(CONVERGENCE_TEST_CTestImplexWrapper)
	, norms(1)
{
}

CTestImplexWrapper::CTestImplexWrapper(ConvergenceTest* test, double maxError,
	double maxReduction, FloorPolicy floor, int print)
	: ConvergenceTest(CONVERGENCE_TEST_CTestImplexWrapper)
	, theTest(test)
	, maxImplexError(maxError)
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
	ConvergenceTest* innerCopy = theTest ? theTest->getCopy(iterations) : 0;
	return new CTestImplexWrapper(innerCopy, maxImplexError, maxReductionFactor,
		onFloor, printFlag);
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

double CTestImplexWrapper::aggregateImplexError(double& minTimeRatio)
{
	// measure this rank's material points. Only the ones that took part in
	// this step, and only those not measured yet: a second call in the same
	// step costs nothing
	const IMPLEXManager::Aggregate& agg = IMPLEXManager::instance().aggregate();
	double e = agg.any_nan ? std::numeric_limits<double>::quiet_NaN() : agg.max;
	minTimeRatio = agg.min_time_ratio;

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
	// THE REDUCTION LIVES HERE, not in the registry, because this is the only
	// place that knows the call is collective: test() is reached by every rank
	// the same number of times, while the registry is just an accumulator with
	// no business knowing about MPI.
	//
	// The guard is the one MPCORecorder uses for the same situation: before
	// the model is partitioned the object exists on P0 only, and a collective
	// would hang. Once partitioned, every rank has its own copy of this test -
	// specifyCTest sends it to the subdomains - so all of them get here.
	//
	// NOTE for OpenSeesSP: this assumes master and subdomain analyses reach
	// test() in lockstep. If they do not, this hangs - loudly, and never with
	// a silently wrong answer, which is why it is worth trying rather than
	// refusing to run. A partial maximum would be the bad outcome: ranks
	// rejecting different steps.
	int pid = 0, np = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	if (np > 1 && !(pid == 0 && !OPS_PARTITIONED)) {
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
	}
#endif

	return e;
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
	double minTimeRatio = 1.0;
	double e = aggregateImplexError(minTimeRatio);
	lastImplexError = e;

	if (std::isnan(e)) {
		// no metric is not a small metric: a material that could not solve
		// implicitly cannot say anything about its own error
		if (printFlag != 0)
			opserr << "CTestImplexWrapper: no IMPL-EX metric (NaN), step rejected\n";
		return -2;
	}

	if (e <= maxImplexError)
		return result;                   // both criteria satisfied

	if (onFloor == Floor_Accept && minTimeRatio <= maxReductionFactor) {
		// the step cannot be reduced any further: taking it is better than
		// deadlocking. The other policy, which is the one to use under an
		// external adaptive scheme, is to fail and let it decide
		if (printFlag != 0)
			opserr << "CTestImplexWrapper: IMPL-EX error " << e << " > " << maxImplexError
			<< " but the step is at the floor (" << minTimeRatio << "), accepted\n";
		return result;
	}

	if (printFlag != 0)
		opserr << "CTestImplexWrapper: IMPL-EX error " << e << " > " << maxImplexError
		<< ", step rejected\n";

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
	// other side can build one to receive into
	static Vector x(5);
	x(0) = maxImplexError;
	x(1) = maxReductionFactor;
	x(2) = static_cast<double>(static_cast<int>(onFloor));
	x(3) = static_cast<double>(printFlag);
	x(4) = static_cast<double>(theTest ? theTest->getClassTag() : -1);
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
	static Vector x(5);
	if (theChannel.recvVector(this->getDbTag(), cTag, x) < 0) {
		opserr << "CTestImplexWrapper::recvSelf() - failed to receive data\n";
		return -1;
	}
	maxImplexError = x(0);
	maxReductionFactor = x(1);
	onFloor = static_cast<FloorPolicy>(static_cast<int>(x(2)));
	printFlag = static_cast<int>(x(3));
	int innerClassTag = static_cast<int>(x(4));
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
