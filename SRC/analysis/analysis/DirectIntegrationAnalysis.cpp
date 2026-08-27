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
                                                                        
// $Revision: 1.15 $
// $Date: 2009-05-14 22:50:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/DirectIntegrationAnalysis.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the implementation of the
// DirectIntegrationAnalysis class.
//
// What: "@(#) DirectIntegrationAnalysis.C, revA"


#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <FE_EleIter.h>

#include <DirectIntegrationAnalysis.h>
#include <EquiSolnAlgo.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <ConvergenceTest.h>
#include <TransientIntegrator.h>
#include <Domain.h>

#include <FE_Element.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>
#include <Matrix.h>
#include <ID.h>
#include <Graph.h>

// Constructor
//    sets theModel and theSysOFEqn to 0 and the Algorithm to the one supplied

DirectIntegrationAnalysis::DirectIntegrationAnalysis(Domain &the_Domain,
						     ConstraintHandler &theHandler,
						     DOF_Numberer &theNumberer,
						     AnalysisModel &theModel,
						     EquiSolnAlgo &theSolnAlgo,		   
						     LinearSOE &theLinSOE,
						     TransientIntegrator &theTransientIntegrator,
						     ConvergenceTest *theConvergenceTest,
						     int num_SubLevels, 
						     int num_SubSteps)
:TransientAnalysis(the_Domain), 
 theConstraintHandler(&theHandler),
 theDOF_Numberer(&theNumberer), 
 theAnalysisModel(&theModel), 
 theAlgorithm(&theSolnAlgo), 
 theSOE(&theLinSOE), 
 theEigenSOE(0),
 theIntegrator(&theTransientIntegrator), 
 theTest(theConvergenceTest),
 domainStamp(0),
 numSubLevels(num_SubLevels),
 numSubSteps(num_SubSteps)
{
  // first we set up the links needed by the elements in the 
  // aggregation
  theAnalysisModel->setLinks(the_Domain, theHandler);
  theConstraintHandler->setLinks(the_Domain, theModel, theTransientIntegrator);
  theDOF_Numberer->setLinks(theModel);
  theIntegrator->setLinks(theModel, theLinSOE, theTest);
  theAlgorithm->setLinks(theModel, theTransientIntegrator, theLinSOE, theTest);
  theSOE->setLinks(theModel);

  if (theTest != 0)
    theAlgorithm->setConvergenceTest(theTest);
  else
    theTest = theAlgorithm->getConvergenceTest();
  
}    

DirectIntegrationAnalysis::~DirectIntegrationAnalysis()
{
  // we don't invoke the destructors in case user switching
  // from a static to a direct integration analysis 
  // clearAll() must be invoked if user wishes to invoke destructor
}    

void
DirectIntegrationAnalysis::clearAll(void)
{
  // invoke the destructor on all the objects in the aggregation
  if (theAnalysisModel != 0)     
    delete theAnalysisModel;
  if (theConstraintHandler != 0) 
    delete theConstraintHandler;
  if (theDOF_Numberer != 0)      
    delete theDOF_Numberer;
  if (theIntegrator != 0) 
    delete theIntegrator;
  if (theAlgorithm != 0)  
    delete theAlgorithm;
  if (theSOE != 0)
    delete theSOE;
  if (theEigenSOE != 0)
    delete theEigenSOE;
  if (theTest != 0)
    delete theTest;


    theAnalysisModel =0;
    theConstraintHandler =0;
    theDOF_Numberer =0;
    theIntegrator =0;
    theAlgorithm =0;
    theSOE =0;
    theEigenSOE =0;
    theTest =0;
}    

#include <NodeIter.h>
#include <Node.h>

int 
DirectIntegrationAnalysis::initialize(void)
{
    Domain *the_Domain = this->getDomainPtr();

    // check if domain has undergone change
    int stamp = the_Domain->hasDomainChanged();
    if (this->anyDomainChange(stamp != domainStamp)) {
      domainStamp = stamp;	
      if (this->domainChanged() < 0) {
	opserr << "DirectIntegrationAnalysis::initialize() - domainChanged() failed\n";
	return -1;
      }	
    }
    if (theIntegrator->initialize() < 0) {
	opserr << "DirectIntegrationAnalysis::initialize() - integrator initialize() failed\n";
	return -2;
    } else
      theIntegrator->commit();

    return 0;
}


int 
DirectIntegrationAnalysis::analyze(int numSteps, double dT, bool flush)
{
  int result = 0;

  for (int i=0; i<numSteps; i++) {
    result = this->analyzeStep(dT);
    if (result < 0) {
      if (numSubLevels != 0)
	result = this->analyzeSubLevel(1, dT);
      if (result < 0)
	return result;
    }
  }

  Domain *the_Domain = this->getDomainPtr();
  if (the_Domain != 0 && flush) {
    the_Domain->flushRecorders();
  }

  return result;
}

int 
DirectIntegrationAnalysis::analyzeStep(double dT)
{
  int result = 0;
  Domain *the_Domain = this->getDomainPtr();

  // Each phase agrees on its outcome across the processes before anyone acts on
  // it - see Analysis::worstStepResult(). The agreement matters twice over here:
  // besides the collectives of the next step, analyze() and analyzeSubLevel()
  // decide from THIS return value whether to sub-step, and that decision recurses.
  // A local verdict would have one process sub-stepping (more collectives) while
  // another moved on to the next step. The agreement is unconditional.
  // The failure messages below are reported once, by one process, and name the
  // process that failed - see Analysis::reportHere(). Sub-stepping makes a failure
  // here an ordinary event, so N copies of it per rejected step is real noise.
  if (this->worstStepResult(theAnalysisModel->analysisStep(dT), "analysisStep()") < 0) {
    if (this->reportHere()) {
      opserr << "DirectIntegrationAnalysis::analyze() - the AnalysisModel failed";
      opserr << " at time " << the_Domain->getCurrentTime() << this->whoFailed() << endln;
    }
    the_Domain->revertToLastCommit();
    return -2;
  }

  // check if domain has undergone change
  int stamp = the_Domain->hasDomainChanged();
  if (this->anyDomainChange(stamp != domainStamp)) {
    domainStamp = stamp;
    if (this->worstStepResult(this->domainChanged(), "domainChanged()") < 0) {
      if (this->reportHere())
	opserr << "DirectIntegrationAnalysis::analyze() - domainChanged() failed"
	       << this->whoFailed() << endln;
      return -1;
    }
  }

  if (this->worstStepResult(theIntegrator->newStep(dT), "newStep()") < 0) {
    if (this->reportHere()) {
      opserr << "DirectIntegrationAnalysis::analyze() - the Integrator failed";
      opserr << " at time " << the_Domain->getCurrentTime() << this->whoFailed() << endln;
    }
    the_Domain->revertToLastCommit();
    theIntegrator->revertToLastStep();
    return -2;
  }

  result = theAlgorithm->solveCurrentStep();
  result = this->worstStepResult(result, "solveCurrentStep()");
  if (result < 0) {
    if (this->reportHere()) {
      opserr << "DirectIntegrationAnalysis::analyze() - the Algorithm failed";
      opserr << " at time " << the_Domain->getCurrentTime() << this->whoFailed() << endln;
    }
    the_Domain->revertToLastCommit();
    theIntegrator->revertToLastStep();
    return -3;
  }
  
  // AddingSensitivity:BEGIN ////////////////////////////////////
#ifdef _RELIABILITY

    if (theIntegrator->shouldComputeAtEachStep()) {
	
      result = theIntegrator->computeSensitivities();
      if (result < 0) {
	opserr << "DirectIntegrationAnalysis::analyze() - the SensitivityAlgorithm failed";
	opserr << " at time ";
	opserr << the_Domain->getCurrentTime() << endln;
	the_Domain->revertToLastCommit();	    
	theIntegrator->revertToLastStep();
	return -5;
      }    
    }
#endif
  // AddingSensitivity:END //////////////////////////////////////
  
    if (AnalysisCommitFilter::instance().isActive()) {
        result = AnalysisCommitFilter::instance().test();
        // the filter evaluates a Tcl expression, so its verdict is rank-local by
        // nature - the one phase here whose failure nothing else can reduce
        result = this->worstStepResult(result, "commitFilter()");
        if (result < 0) {
            opserr << "DirectIntegrationAnalysis::analyze() - the commit filter failed";
            opserr << " at time " << the_Domain->getCurrentTime() << endln;
            the_Domain->revertToLastCommit();
            theIntegrator->revertToLastStep();
            return -6;
        }
    }

  result = theIntegrator->commit();
  result = this->worstStepResult(result, "commit()");
  if (result < 0) {
    if (this->reportHere()) {
      opserr << "DirectIntegrationAnalysis::analyze() - ";
      opserr << "the Integrator failed to commit";
      opserr << " at time " << the_Domain->getCurrentTime() << this->whoFailed() << endln;
    }
    the_Domain->revertToLastCommit();
    theIntegrator->revertToLastStep();
    return -4;
  } 
    
  return result;
}

int
DirectIntegrationAnalysis::analyzeSubLevel(int level, double dT) {
  int result = 0;
  if (numSubSteps == 0)
    return -1;

  double stepDT = dT/(numSubSteps*1.);

  for (int i=0; i<numSubSteps; i++) {
    result = this->analyzeStep(stepDT);
    if (result < 0) {
      if (level == numSubLevels) {
	return result;
      } else {
	result = this->analyzeSubLevel(level+1, stepDT);
	if (result < 0)
	  return result;
      }
    }
  }
  return result;
}

int 
DirectIntegrationAnalysis::eigen(int numMode, bool generalized, bool findSmallest)
{
    if (theAnalysisModel == 0 || theEigenSOE == 0) {
      opserr << "WARNING DirectIntegrationAnalysis::eigen() - no EigenSOE has been set\n";
      return -1;
    }

    int result = 0;
    Domain *the_Domain = this->getDomainPtr();

    // Every phase below agrees on its outcome across the processes before anyone acts on
    // it - see Analysis::worstStepResult(). theEigenSOE->solve() is a collective in the
    // multi-rank eigen, so a process returning early from an assembly failure leaves
    // the others waiting inside ARPACK with no participant.
    result = this->worstStepResult(theAnalysisModel->eigenAnalysis(numMode, generalized,
								  findSmallest),
				   "eigenAnalysis()");
    if (result < 0) {
      // this return value was assigned and never tested
      if (this->reportHere())
	opserr << "DirectIntegrationAnalysis::eigen() - the AnalysisModel failed"
	       << this->whoFailed() << endln;
      return -1;
    }

    int stamp = the_Domain->hasDomainChanged();

    if (this->anyDomainChange(stamp != domainStamp)) {
      domainStamp = stamp;

      result = this->worstStepResult(this->domainChanged(), "domainChanged()");

      if (result < 0) {
	if (this->reportHere())
	  opserr << "DirectIntegrationAnalysis::eigen() - domainChanged failed"
		 << this->whoFailed() << endln;
	return -1;
      }
    }


    //
    // zero A and M
    //
    theEigenSOE->zeroA();
    theEigenSOE->zeroM();

    //
    // form K
    //


    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    FE_Element *elePtr;

    int numFailedK = 0;
    while((elePtr = theEles()) != 0) {
      elePtr->zeroTangent();
      elePtr->addKtToTang(1.0);
      if (theEigenSOE->addA(elePtr->getTangent(0), elePtr->getID()) < 0) {
	numFailedK++;
	// local detail, so NOT gated: each process names the equations it could not
	// assemble. First in full then a count, as in Domain::update.
	if (numFailedK == 1) {
	  opserr << "WARNING DirectIntegrationAnalysis::eigen() -";
	  opserr << " failed in addA for ID " << elePtr->getID();
	}
	result = -2;
      }
    }
    if (numFailedK > 1)
      opserr << "WARNING DirectIntegrationAnalysis::eigen() - addA failed for "
	     << numFailedK << " elements on this process\n";

    // used to be set and never read before the final `return 0`
    result = this->worstStepResult(result, "eigen formK");
    if (result < 0) {
      if (this->reportHere())
	opserr << "DirectIntegrationAnalysis::eigen() - the stiffness assembly failed"
	       << this->whoFailed() << endln;
      return -2;
    }

    //
    // if generalized is true, form M
    //

    if (generalized == true) {
      // the `int result = 0;` that used to shadow the function's result here was
      // already commented out in this copy; StaticAnalysis::eigen() still had it live
      int numFailedM = 0;
      FE_EleIter &theEles2 = theAnalysisModel->getFEs();
      while((elePtr = theEles2()) != 0) {
	elePtr->zeroTangent();
	elePtr->addMtoTang(1.0);
	if (theEigenSOE->addM(elePtr->getTangent(0), elePtr->getID()) < 0) {
	  numFailedM++;
	  if (numFailedM == 1) {
	    opserr << "WARNING DirectIntegrationAnalysis::eigen() -";
	    // said "addA" here, in the addM loop
	    opserr << " failed in addM for element ID " << elePtr->getID();
	  }
	  result = -2;
	}
      }

      DOF_Group *dofPtr;
      DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();
      while((dofPtr = theDofs()) != 0) {
	dofPtr->zeroTangent();
	dofPtr->addMtoTang(1.0);
	if (theEigenSOE->addM(dofPtr->getTangent(0),dofPtr->getID()) < 0) {
	  numFailedM++;
	  if (numFailedM == 1) {
	    opserr << "WARNING DirectIntegrationAnalysis::eigen() -";
	    opserr << " failed in addM for DOF group ID " << dofPtr->getID();
	  }
	  result = -3;
	}
      }

      if (numFailedM > 1)
	opserr << "WARNING DirectIntegrationAnalysis::eigen() - addM failed "
	       << numFailedM << " times on this process\n";
    }

    result = this->worstStepResult(result, "eigen formM");
    if (result < 0) {
      if (this->reportHere())
	opserr << "DirectIntegrationAnalysis::eigen() - the mass assembly failed"
	       << this->whoFailed() << endln;
      return result;
    }

    //
    // solve for the eigen values & vectors
    //

    result = this->worstStepResult(theEigenSOE->solve(numMode, generalized, findSmallest),
				   "eigenSolve()");
    if (result < 0) {
	if (this->reportHere())
	  opserr << "WARNING DirectIntegrationAnalysis::eigen() - EigenSOE failed in solve()"
		 << this->whoFailed() << endln;
	return -4;
    }
	
    //
    // now set the eigenvalues and eigenvectors in the model
    //

    theAnalysisModel->setNumEigenvectors(numMode);
    Vector theEigenvalues(numMode);
    for (int i = 1; i <= numMode; i++) {
      theEigenvalues[i-1] = theEigenSOE->getEigenvalue(i);
      theAnalysisModel->setEigenvector(i, theEigenSOE->getEigenvector(i));
    }    
    theAnalysisModel->setEigenvalues(theEigenvalues);
  
    return 0;
}



int
DirectIntegrationAnalysis::domainChanged(void)
{
    Domain *the_Domain = this->getDomainPtr();
    int stamp = the_Domain->hasDomainChanged();
    domainStamp = stamp;

    theAnalysisModel->clearAll();    
    theConstraintHandler->clearAll();
    
    // now we invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.
    theConstraintHandler->handle();

    // we now invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.
    theDOF_Numberer->numberDOF();

    theConstraintHandler->doneNumberingDOF();

    // we invoke setGraph() on the LinearSOE which
    // causes that object to determine its size
    Graph &theGraph = theAnalysisModel->getDOFGraph();

    int result = theSOE->setSize(theGraph);
    if (result < 0) {
	opserr << "DirectIntegrationAnalysis::handle() - ";
	opserr << "LinearSOE::setSize() failed";
	return -3;
    }	    

    if (theEigenSOE != 0) {
      result = theEigenSOE->setSize(theGraph);
      if (result < 0) {
	opserr << "DirectIntegrationAnalysis::handle() - ";
	opserr << "EigenSOE::setSize() failed";
	return -3;
      }	    
    }

    theAnalysisModel->clearDOFGraph();

    // we invoke domainChange() on the integrator and algorithm
    theIntegrator->domainChanged();
    theAlgorithm->domainChanged();

    return 0;
}    

int 
DirectIntegrationAnalysis::setNumberer(DOF_Numberer &theNewNumberer) 
{
    // invoke the destructor on the old one
    if (theDOF_Numberer != 0)
	delete theDOF_Numberer;

    // first set the links needed by the Algorithm
    theDOF_Numberer = &theNewNumberer;
    theDOF_Numberer->setLinks(*theAnalysisModel);

    // invoke domainChanged() either indirectly or directly
    domainStamp = 0;
    return 0;
}



int 
DirectIntegrationAnalysis::setAlgorithm(EquiSolnAlgo &theNewAlgorithm) 
{
  // invoke the destructor on the old one
  if (theAlgorithm != 0)
    delete theAlgorithm;
  
  // first set the links needed by the Algorithm
  theAlgorithm = &theNewAlgorithm;

  if (theAnalysisModel != 0 && theIntegrator != 0 && theSOE != 0)
    theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);
  // invoke domainChanged() either indirectly or directly
  // domainStamp = 0;
  if (domainStamp != 0)
    theAlgorithm->domainChanged();

  return 0;
}


int 
DirectIntegrationAnalysis::setIntegrator(TransientIntegrator &theNewIntegrator)
{
  // invoke the destructor on the old one
  if (theIntegrator != 0) {
      delete theIntegrator;
  }
  // set the links needed by the other objects in the aggregation
  Domain *the_Domain = this->getDomainPtr();
  theIntegrator = &theNewIntegrator;
  theIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);
  theConstraintHandler->setLinks(*the_Domain, *theAnalysisModel, *theIntegrator);
  theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);

  // cause domainChanged to be invoked on next analyze
  //  domainStamp = 0;
  if (domainStamp != 0)
    theIntegrator->domainChanged();
   
  return 0;
}
int 
DirectIntegrationAnalysis::setLinearSOE(LinearSOE &theNewSOE)
{
  // invoke the destructor on the old one
  if (theSOE != 0)
    delete theSOE;

  // set the links needed by the other objects in the aggregation
  theSOE = &theNewSOE;
  theIntegrator->setLinks(*theAnalysisModel,*theSOE, theTest);
  theAlgorithm->setLinks(*theAnalysisModel, *theIntegrator, *theSOE, theTest);
  theSOE->setLinks(*theAnalysisModel);
  
  if (theEigenSOE != 0) 
    theEigenSOE->setLinearSOE(*theSOE);
  
  // cause domainChanged to be invoked on next analyze
  domainStamp = 0;
  
  return 0;
}

int 
DirectIntegrationAnalysis::setEigenSOE(EigenSOE &theNewSOE)
{
  // invoke the destructor on the old one if not the same!
  if (theEigenSOE != 0) {
    if (theEigenSOE->getClassTag() != theNewSOE.getClassTag()) {
      delete theEigenSOE;
      theEigenSOE = 0;
    }
  }

  if (theEigenSOE == 0) {
    theEigenSOE = &theNewSOE;
    theEigenSOE->setLinks(*theAnalysisModel);
    theEigenSOE->setLinearSOE(*theSOE);

    /*
    if (domainStamp != 0) {
      Graph &theGraph = theAnalysisModel->getDOFGraph();
      theEigenSOE->setSize(theGraph);
    }
    */
    domainStamp = 0;
  }
 
  return 0;
}

int 
DirectIntegrationAnalysis::setConvergenceTest(ConvergenceTest &theNewTest)
{
  // invoke the destructor on the old one
  if (theTest != 0)
    delete theTest;
  
  // set the links needed by the other objects in the aggregation
  theTest = &theNewTest;

  if (theIntegrator != 0)
    theIntegrator->setLinks(*theAnalysisModel, *theSOE, theTest);

  if (theAlgorithm != 0)
    theAlgorithm->setConvergenceTest(theTest);
  
  return 0;
}


int
DirectIntegrationAnalysis::checkDomainChange(void)
{
  Domain *the_Domain = this->getDomainPtr();

  // check if domain has undergone change
  int stamp = the_Domain->hasDomainChanged();
  if (this->anyDomainChange(stamp != domainStamp)) {
    domainStamp = stamp;
    if (this->domainChanged() < 0) {
      opserr << "DirectIntegrationAnalysis::initialize() - domainChanged() failed\n";
      return -1;
    }
  }

  return 0;
}


EquiSolnAlgo *
DirectIntegrationAnalysis::getAlgorithm(void)
{
  return theAlgorithm;
}

AnalysisModel *
DirectIntegrationAnalysis::getModel(void)
{
  return theAnalysisModel;
}


TransientIntegrator *
DirectIntegrationAnalysis::getIntegrator(void)
{
  return theIntegrator;
}

ConvergenceTest *
DirectIntegrationAnalysis::getConvergenceTest(void)
{
  return theTest;
}




