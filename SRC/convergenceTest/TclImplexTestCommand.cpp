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
// The Tcl 'implexTest' command: wraps the convergence test that is currently
// installed, so that the IMPL-EX error criterion is ANDed on top of it.
//
//     test NormDispIncr 1.0e-6 10 0
//     implexTest -tol 0.05
//
//     implexTest <-tol $tol> <-maxReduction $f> <-onFloor fail|accept>
//                <-print $flag>
//
// WHY A COMMAND OF ITS OWN, and not another type of 'test'. 'test' is
// positional after the type name and builds one object out of its own
// arguments; this one builds nothing of its own, it decorates what is already
// there. Keeping the two apart means specifyCTest is not touched at all - only
// two lines of commands.cpp are, to register this command - which is what makes
// the feature survive an upstream rewrite of that function without a conflict.
//
// WHY IT COPIES THE TEST IT WRAPS. The test being decorated is already
// installed, and installing the wrapper deletes whatever it replaces - which is
// exactly that test. So the wrapper is handed a getCopy() of it and the original
// is left for the install to dispose of.
//
// IT MUST FOLLOW EVERY 'test'. Declaring a test replaces the wrapper, and there
// is no way for this command to notice and reattach itself. Whoever writes the
// input is responsible for emitting the two together.

#include <tcl.h>
#include <string.h>

#include <OPS_Globals.h>
#include <ConvergenceTest.h>
#include <CTestImplexWrapper.h>
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <classTags.h>

#ifdef _PARALLEL_PROCESSING
#include <Domain.h>
#include <Subdomain.h>
#include <SubdomainIter.h>
extern Domain theDomain;
#endif

// the analysis state of commands.cpp. This command replaces the installed
// convergence test in exactly the way 'test' does, so it needs the same three
extern ConvergenceTest *theTest;
extern StaticAnalysis *theStaticAnalysis;
extern DirectIntegrationAnalysis *theTransientAnalysis;

int
TclImplexTestCommand(ClientData clientData, Tcl_Interp *interp, int argc,
		     TCL_Char **argv)
{
  if (theTest == 0) {
    opserr << "WARNING implexTest - there is no convergence test to wrap: "
	   << "declare one with 'test' first\n";
    return TCL_ERROR;
  }

  // wrapping a wrapper would give two IMPL-EX criteria in series, each clearing
  // the registry the other is about to read
  if (theTest->getClassTag() == CONVERGENCE_TEST_CTestImplexWrapper) {
    opserr << "WARNING implexTest - the current convergence test is already an "
	   << "IMPL-EX wrapper. Re-declare the inner test with 'test' before "
	   << "calling this again\n";
    return TCL_ERROR;
  }

  double implexTol = 0.05;
  double maxReduction = 0.0;
  int floorPolicy = CTestImplexWrapper::Floor_Fail;
  int printImplex = 0;

  int i = 1;
  while (i < argc) {

    if (strcmp(argv[i],"-tol") == 0) {
      if (i + 1 >= argc) {
	opserr << "WARNING implexTest - -tol needs a value\n";
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[++i], &implexTol) != TCL_OK)
	return TCL_ERROR;
      if (implexTol <= 0.0) {
	opserr << "WARNING implexTest - -tol must be positive, got "
	       << implexTol << "\n";
	return TCL_ERROR;
      }
      ++i;

    } else if (strcmp(argv[i],"-maxReduction") == 0) {
      if (i + 1 >= argc) {
	opserr << "WARNING implexTest - -maxReduction needs a value\n";
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[++i], &maxReduction) != TCL_OK)
	return TCL_ERROR;
      if (maxReduction < 0.0 || maxReduction >= 1.0) {
	opserr << "WARNING implexTest - -maxReduction is a fraction of the "
	       << "nominal time step and must be in [0,1), got "
	       << maxReduction << "\n";
	return TCL_ERROR;
      }
      ++i;

    } else if (strcmp(argv[i],"-onFloor") == 0) {
      if (i + 1 >= argc) {
	opserr << "WARNING implexTest - -onFloor needs fail or accept\n";
	return TCL_ERROR;
      }
      ++i;
      if (strcmp(argv[i],"accept") == 0)
	floorPolicy = CTestImplexWrapper::Floor_Accept;
      else if (strcmp(argv[i],"fail") == 0)
	floorPolicy = CTestImplexWrapper::Floor_Fail;
      else {
	opserr << "WARNING implexTest - -onFloor takes fail or accept, got "
	       << argv[i] << "\n";
	return TCL_ERROR;
      }
      ++i;

    } else if (strcmp(argv[i],"-print") == 0) {
      if (i + 1 >= argc) {
	opserr << "WARNING implexTest - -print needs a value\n";
	return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[++i], &printImplex) != TCL_OK)
	return TCL_ERROR;
      ++i;

    } else {
      opserr << "WARNING implexTest - unknown option " << argv[i] << "\n";
      return TCL_ERROR;
    }
  }

  // the wrapper takes ownership of this copy. The original cannot be handed
  // over directly: it is the test the install below is about to delete
  ConvergenceTest *theInnerTest = theTest->getCopy(theTest->getMaxNumTests());
  if (theInnerTest == 0) {
    opserr << "WARNING implexTest - the current convergence test could not be "
	   << "copied\n";
    return TCL_ERROR;
  }

  theTest = new CTestImplexWrapper(theInnerTest, implexTol, maxReduction,
      static_cast<CTestImplexWrapper::FloorPolicy>(floorPolicy), printImplex);

  // from here on this is the install of specifyCTest, verbatim. Note what it
  // does NOT do: when no analysis exists yet the test being replaced is not
  // deleted, because nobody owns it and it cannot be deleted here either - an
  // algorithm declared before this point was given it by reference. A second
  // 'test' leaves it behind in the same way

  // if the analysis exists - we want to change the Test
  if (theStaticAnalysis != 0)
    theStaticAnalysis->setConvergenceTest(*theTest);

  else if (theTransientAnalysis != 0)
    theTransientAnalysis->setConvergenceTest(*theTest);

#ifdef _PARALLEL_PROCESSING
  if (theStaticAnalysis != 0 || theTransientAnalysis != 0) {
    SubdomainIter &theSubdomains = theDomain.getSubdomains();
    Subdomain *theSub;
    while ((theSub = theSubdomains()) != 0) {
      theSub->setAnalysisConvergenceTest(*theTest);
    }
  }
#endif

  return TCL_OK;
}
