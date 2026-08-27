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
// A convergence test that decorates another one and ANDs the IMPL-EX error
// criterion on top of it. Built by the 'implexTest' command, which wraps
// whatever test is installed at the time:
//
//     test NormDispIncr 1.0e-6 10 0
//     implexTest -tol 0.05 <-maxFraction 0.0> <-maxReduction 0.01>
//                <-onFloor fail|accept> <-print 0>
//
// See TclImplexTestCommand.cpp - in particular for why it must follow EVERY
// 'test', and why it wraps a copy.
//
// WHAT THE CRITERION IS. Not 'the worst point is within tolerance' but 'few
// enough points are over it':
//
//     reject  <=>  #(e_i > tol) / N  >  maxFraction
//
// with N the points that took part in the step - not the points of the model,
// so an element deactivated by staged construction is not in the denominator.
// maxFraction defaults to ZERO, and at zero this IS the worst-point criterion:
// the count is nonzero exactly when the maximum is over tolerance. Above zero
// it stops a single gauss point out of millions from driving the step size of
// the whole model down by itself.
//
// A fraction and not a mean, because a mean over material points answers a
// different question - see the note on count_over in IMPLEXManager.h. And a
// fraction has one consequence worth stating: with maxFraction > 0 an ISOLATED
// error of any size passes, because it is isolated. There is deliberately no
// hard cap on the maximum; the maximum stays readable through
// getImplexError(), so a scheme that wants one can put it outside.
//
// WHY A WRAPPER. An analysis has other nonlinearities than the extrapolation,
// and the test the user chose is what knows about them. So that one runs
// first and can already say 'not converged'; only when it says converged does
// the IMPL-EX criterion get a say, and the answer is the AND of the two.
//
// WHY THE MATERIAL DOES NOT DO THIS. The error is local, the decision is not:
// a material that fails the step owns a policy it cannot see the analysis to
// choose, and the code it returns reaches the element as a material failure,
// not as a controlled rejection. The materials measure; this decides.
//
// WHEN IT MEASURES: only on the iteration where the inner test says
// converged. Before that the state is not equilibrated and the metric means
// nothing, and measuring at every iteration would pay one implicit solve per
// iteration per material point, which is exactly what this design exists to
// avoid. Every equiSolnAlgo calls test() inside a 'while (result == -1)'
// loop, so this happens ONCE per solveCurrentStep - once per attempted step,
// not per iteration.
//
// WHY IT RETURNS -2 AND NOT -1 ON A REJECTION. -1 means 'keep iterating', and
// iterating cannot fix this: the strain is what it is and the error is one of
// discretisation, so Newton would grind to maxNumIter and fail anyway, having
// paid for it. -2 is 'failed', which is what StaticAnalysis turns into a
// revertToLastCommit + revertToLastStep - the step is undone, and the adaptive
// scheme outside (in Tcl) halves the step and tries again.
//
// WHAT IT DOES NOT DO: propose a new step size. A convergence test cannot
// reduce the step and does not need to. Measured, on the same material,
// protocol and tolerances: a fixed halving with growth capped at the nominal
// step is not worse than the proportional factor safety*(tol/e)^(1/q) - same
// accuracy within 5% with 10-19% FEWER steps across the whole tolerance
// sweep - because the error here is bimodal and concentrated at the events
// (exactly zero on 511 of 600 steps), while the proportional model assumes it
// is smooth. So this is a yes/no gate, and the schedule stays in Tcl.

#ifndef CTestImplexWrapper_h
#define CTestImplexWrapper_h

#include <ConvergenceTest.h>
#include <Vector.h>

class EquiSolnAlgo;

class CTestImplexWrapper : public ConvergenceTest
{
public:
	// what to do when the step cannot be reduced any further
	enum FloorPolicy {
		// fail anyway, and let whoever drives the analysis decide. This is the
		// one to use with an external adaptive scheme, which has its own
		// maximum reduction and its own idea of what to do at it
		Floor_Fail = 0,
		// accept the step even though it is over tolerance, rather than
		// deadlock
		Floor_Accept = 1
	};

public:
	// constructors
	CTestImplexWrapper();
	// takes OWNERSHIP of theTest
	CTestImplexWrapper(ConvergenceTest* theTest, double maxImplexError,
		double maxFraction, double maxReductionFactor, FloorPolicy onFloor,
		int printFlag);

	// destructor
	~CTestImplexWrapper();

	ConvergenceTest* getCopy(int iterations);

	int setEquiSolnAlgo(EquiSolnAlgo& theAlgo);

	int test(void);
	int start(void);

	int getNumTests(void);
	int getMaxNumTests(void);
	double getRatioNumToMax(void);
	// the norms of the wrapped test with the IMPL-EX error appended as the
	// LAST component. It is the only channel a Tcl script has to read the
	// error after a failure, which is what an adaptive scheme wants
	const Vector& getNorms(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// the largest IMPL-EX error of the last measured step. Survives the
	// failure it caused, on purpose
	inline double getImplexError(void) const { return lastImplexError; }

private:
	// what one measurement says, AFTER the reduction over the ranks
	struct Measure {
		// the largest error anywhere in the model, or NaN if some point could
		// not be measured at all
		double max = 0.0;
		// of the points that took part, the fraction that is over tolerance.
		// Zero when none is, including when none took part
		double fraction = 0.0;
		// the two counts behind it. Only the diagnostic reads them, and it
		// must: 'over 1 of 240000' and 'over 900 of 240000' are the same
		// fraction to nobody
		double over = 0.0;
		double count = 0.0;
		// how far the schedule has cut, relative to this analysis step's
		// nominal - see nominalDt
		double minTimeRatio = 1.0;
	};

	// measure this rank's points and reduce over all ranks. See the note in the
	// implementation on why the reduction lives here and not in the registry
	Measure aggregateImplexError(void);

private:
	// the test chosen by the user. Owned
	ConvergenceTest* theTest = 0;
	// the largest IMPL-EX error a point may carry before it counts as a
	// violation
	double maxImplexError = 0.05;
	// the fraction of the points that may carry one before the step is
	// rejected. ZERO IS THE DEFAULT AND IS NOT 'disabled': at zero the
	// criterion is 'no point may be over tolerance', which is the worst-point
	// criterion this test has always applied
	double maxFraction = 0.0;
	// the fraction of THIS analysis step's nominal dt under which the step
	// cannot be reduced any further - see nominalDt below
	double maxReductionFactor = 0.0;
	// what to do there
	FloorPolicy onFloor = Floor_Fail;
	// print flag
	int printFlag = 0;
	// the last measured error, kept after the step it failed
	double lastImplexError = 0.0;
	// THE NOMINAL STEP OF THIS ANALYSIS STEP, against which maxReductionFactor
	// is a fraction: the largest pseudo-time increment this object has seen.
	//
	// It lives here, and not in the materials, because it is part of the
	// POLICY - the same reason the tolerance is here. A material's own
	// reference (dtime_0) is frozen at its FIRST commit ever and then outlives
	// every analysis step: a gravity step at dt = 0.1 followed by a
	// displacement-control step at dt = 1.0e-3 has the second one reading
	// ratio = 0.01, so with maxReduction = 0.01 it starts AT the floor and
	// -onFloor accept takes every step of it without ever reducing. Measured
	// on exactly that pair.
	//
	// A running maximum rather than 'the first one seen': the first step of a
	// step may itself be a retry the schedule has already cut, and the nominal
	// is what the schedule grows BACK to. Not serialized - in OpenSeesSP the
	// wrapper travels to the subdomains before any step exists, so both sides
	// start at zero and fill it in from their own first measurement.
	double nominalDt = 0.0;
	// the norms of the inner test plus the error
	Vector norms;
};

#endif
