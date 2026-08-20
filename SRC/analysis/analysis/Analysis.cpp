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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/Analysis.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/analysis/Analysis.C
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the implementation of Analysis.
// Analysis is an abstract base class, i.e. no objects of it's
// type can be created. 
//
// What: "@(#) Analysis.C, revA"

#include <Analysis.h>
#include <Domain.h>
#include <elementAPI.h>

int OPS_SetAnalysisCommitFilter()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if (numData > 0) {
        AnalysisCommitFilter::instance().setExpression(OPS_GetString());
    }
    else {
        AnalysisCommitFilter::instance().unset();
    }
    return 0;
}

AnalysisCommitFilter& AnalysisCommitFilter::instance()
{
    static AnalysisCommitFilter _instance;
    return _instance;
}

void AnalysisCommitFilter::setExpression(const std::string& x)
{
    m_expression = x;
    m_active = true;
}

void AnalysisCommitFilter::unset()
{
    m_active = false;
}

void AnalysisCommitFilter::setCustomFunction(function_t the_custom_function)
{
    m_function = the_custom_function;
}

int AnalysisCommitFilter::test()
{
    return m_function(m_expression);
}

AnalysisCommitFilter::function_t AnalysisCommitFilter::makeDefaultTclFunction()
{
    return [](const std::string& x) -> int {
        double value = 0.0;
        if (OPS_EvalDoubleStringExpression(x.data(), value) < 0)
            return 0;
        return static_cast<int>(value);
    };
}


#include <stdio.h>

Analysis::Analysis(Domain &theDom)
:theDomain(&theDom)
{
  failInfo[0] = '\0';
}

Analysis::~Analysis()
{

}

// for parallel
#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#include <OPS_Globals.h>
#endif // _PARALLEL_INTERPRETERS

int
Analysis::worstStepResult(int resultHere, const char *phase)
{
    failInfo[0] = '\0';

#ifdef _PARALLEL_INTERPRETERS
    int np = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (np > 1) {

      int here = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &here);

      // MINLOC over the pair (code, process) gets the worst code AND its owner out
      // of the one reduce that the agreement needs anyway. Where several processes
      // report the same worst code MPI_MINLOC returns the lowest-numbered of them,
      // so every process names the same process and the diagnostic is the same on
      // every run - unlike whichever stderr line happens to arrive first.
      int in[2], out[2];
      in[0] = resultHere;  in[1] = here;
      out[0] = resultHere; out[1] = here;

      if (MPI_Allreduce(in, out, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD) != MPI_SUCCESS) {
	opserr << "Analysis::worstStepResult - MPI_Allreduce failed\n";
	return resultHere;
      }

      if (out[0] < 0) {

	// Whether one process failed or all of them separates a local problem -
	// one partition's element or material - from a global one, typically the
	// solver. It costs a second reduce, on the failure path only, where the
	// analysis is about to stop or sub-step anyway.
	//
	// Guarding a collective on the reduced value is safe, and only on the
	// REDUCED value: out[0] is identical on every process, so this branch is
	// taken by all of them or by none. Guarding on resultHere would not be.
	int failedHere = (resultHere < 0) ? 1 : 0;
	int numFailed = failedHere;
	if (MPI_Allreduce(&failedHere, &numFailed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
	  numFailed = 0;

	sprintf(failInfo, " [%.24s returned %d on process %d of %d; %d process(es) failed]",
		(phase != 0) ? phase : "the phase", out[0], out[1], np, numFailed);
      }

      return out[0];
    }
#endif // _PARALLEL_INTERPRETERS

    return resultHere;
}

bool
Analysis::reportHere(void) const
{
#ifdef _PARALLEL_INTERPRETERS
    int np = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (np > 1) {
      int here = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &here);
      return (here == 0);
    }
#endif // _PARALLEL_INTERPRETERS

    return true;
}

const char *
Analysis::whoFailed(void) const
{
    return failInfo;
}

bool
Analysis::anyDomainChange(bool domainChangedHere)
{
#ifdef _PARALLEL_INTERPRETERS
    int np = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (np > 1) {
      int here = domainChangedHere ? 1 : 0;
      int anywhere = here;
      if (MPI_Allreduce(&here, &anywhere, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
	opserr << "Analysis::anyDomainChange - MPI_Allreduce failed\n";
	return domainChangedHere;
      }
      return (anywhere != 0);
    }
#endif // _PARALLEL_INTERPRETERS

    return domainChangedHere;
}

Domain *
Analysis::getDomainPtr(void)
{
    return theDomain;
}
