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
#include <ParallelAgreement.h>

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

int
Analysis::worstStepResult(int resultHere, const char *phase)
{
    failInfo[0] = '\0';

    // the agreement, the attribution and the lockstep check are all in
    // OPS_agreeWorstWho() - see ParallelAgreement.h. All this adds is the wording,
    // and the out parameters keep their values when the run is not coupled.
    int worstRank = 0, numFailed = 0, np = 1;
    int worst = OPS_agreeWorstWho(resultHere, worstRank, numFailed, np);

    if (worst < 0 && np > 1)
      sprintf(failInfo, " [%.24s returned %d on process %d of %d; %d process(es) failed]",
	      (phase != 0) ? phase : "the phase", worst, worstRank, np, numFailed);

    return worst;
}

bool
Analysis::reportHere(void) const
{
    return (OPS_agreementRank() == 0);
}

const char *
Analysis::whoFailed(void) const
{
    return failInfo;
}

bool
Analysis::anyDomainChange(bool domainChangedHere)
{
    return OPS_agreeAny(domainChangedHere);
}

Domain *
Analysis::getDomainPtr(void)
{
    return theDomain;
}
