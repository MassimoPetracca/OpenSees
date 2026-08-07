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
// The interpreter side of 'implexTest': wraps the convergence test that is
// currently installed, so that the IMPL-EX error criterion is ANDed on top of
// it.
//
//     ops.test('NormDispIncr', 1.0e-6, 10, 0)
//     ops.implexTest('-tol', 0.05)
//
//     implexTest <-tol tol> <-maxReduction f> <-onFloor 'fail'|'accept'>
//                <-print flag>
//
// See TclImplexTestCommand.cpp for why this is a command of its own and why it
// copies the test it wraps. One difference from the Tcl side, and it is in this
// one's favour: OpenSeesCommands::setCTest() deletes the test it replaces
// whether or not an analysis exists, so nothing is left behind here.

#include <string.h>

#include <OpenSeesCommands.h>
#include <ConvergenceTest.h>
#include <CTestImplexWrapper.h>
#include <classTags.h>
#include <elementAPI.h>

int OPS_ImplexCTest()
{
    OpenSeesCommands* cmds = OPS_GetOpenSeesCommands();
    if (cmds == 0) {
	opserr << "WARNING implexTest - no model\n";
	return -1;
    }

    ConvergenceTest* theCurrentTest = cmds->getCTest();
    if (theCurrentTest == 0) {
	opserr << "WARNING implexTest - there is no convergence test to wrap: "
	       << "declare one with test() first\n";
	return -1;
    }

    // wrapping a wrapper would give two IMPL-EX criteria in series, each
    // clearing the registry the other is about to read
    if (theCurrentTest->getClassTag() == CONVERGENCE_TEST_CTestImplexWrapper) {
	opserr << "WARNING implexTest - the current convergence test is already "
	       << "an IMPL-EX wrapper. Re-declare the inner test with test() "
	       << "before calling this again\n";
	return -1;
    }

    double implexTol = 0.05;
    double maxReduction = 0.0;
    int floorPolicy = CTestImplexWrapper::Floor_Fail;
    int printImplex = 0;
    int numData = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {

	const char* opt = OPS_GetString();

	if (strcmp(opt, "-tol") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING implexTest - -tol needs a value\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&numData, &implexTol) < 0) {
		opserr << "WARNING implexTest - failed to read -tol\n";
		return -1;
	    }
	    if (implexTol <= 0.0) {
		opserr << "WARNING implexTest - -tol must be positive, got "
		       << implexTol << "\n";
		return -1;
	    }

	} else if (strcmp(opt, "-maxReduction") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING implexTest - -maxReduction needs a value\n";
		return -1;
	    }
	    if (OPS_GetDoubleInput(&numData, &maxReduction) < 0) {
		opserr << "WARNING implexTest - failed to read -maxReduction\n";
		return -1;
	    }
	    if (maxReduction < 0.0 || maxReduction >= 1.0) {
		opserr << "WARNING implexTest - -maxReduction is a fraction of "
		       << "the nominal time step and must be in [0,1), got "
		       << maxReduction << "\n";
		return -1;
	    }

	} else if (strcmp(opt, "-onFloor") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING implexTest - -onFloor needs fail or accept\n";
		return -1;
	    }
	    const char* value = OPS_GetString();
	    if (strcmp(value, "accept") == 0)
		floorPolicy = CTestImplexWrapper::Floor_Accept;
	    else if (strcmp(value, "fail") == 0)
		floorPolicy = CTestImplexWrapper::Floor_Fail;
	    else {
		opserr << "WARNING implexTest - -onFloor takes fail or accept, "
		       << "got " << value << "\n";
		return -1;
	    }

	} else if (strcmp(opt, "-print") == 0) {
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING implexTest - -print needs a value\n";
		return -1;
	    }
	    if (OPS_GetIntInput(&numData, &printImplex) < 0) {
		opserr << "WARNING implexTest - failed to read -print\n";
		return -1;
	    }

	} else {
	    opserr << "WARNING implexTest - unknown option " << opt << "\n";
	    return -1;
	}
    }

    // the wrapper takes ownership of this copy. The original cannot be handed
    // over directly: it is the test setCTest is about to delete
    ConvergenceTest* theInnerTest =
	theCurrentTest->getCopy(theCurrentTest->getMaxNumTests());
    if (theInnerTest == 0) {
	opserr << "WARNING implexTest - the current convergence test could not "
	       << "be copied\n";
	return -1;
    }

    cmds->setCTest(new CTestImplexWrapper(theInnerTest, implexTol, maxReduction,
	static_cast<CTestImplexWrapper::FloorPolicy>(floorPolicy), printImplex));

    return 0;
}
