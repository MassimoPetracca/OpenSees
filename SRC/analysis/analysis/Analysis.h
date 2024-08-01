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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/Analysis.h,v $
                                                                        
                                                                        
#ifndef Analysis_h
#define Analysis_h

// File: ~/analysis/analysis/Analysis.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the interface for the Analysis class.
// Analysis is an abstract class, i.e. no objects of it's type can be created. 
//
// What: "@(#) Analysis.h, revA"


#include <string>
#include <functional>
class AnalysisCommitFilter
{
public:
    using function_t = std::function<int(const std::string&)>;
private:
    AnalysisCommitFilter() = default;
    AnalysisCommitFilter(const AnalysisCommitFilter&) = delete;
    AnalysisCommitFilter& operator = (const AnalysisCommitFilter&) = delete;
public:
    static AnalysisCommitFilter& instance();
    void setExpression(const std::string& x);
    void unset();
    void setCustomFunction(function_t the_custom_function);
    int test();
    inline bool isActive() const { return m_active; }
private:
    static function_t makeDefaultTclFunction();
private:
    bool m_active = false;
    std::string m_expression;
    function_t m_function = makeDefaultTclFunction();
};

class Domain;

class Analysis
{
  public:
    Analysis(Domain &theDomain);
    virtual ~Analysis();

    // pure virtual functions
    //    virtual int analyze(void) =0;
    virtual int domainChanged(void) = 0;
    
  protected:
    Domain *getDomainPtr(void);
    
  private:
    Domain *theDomain;
};

#endif

