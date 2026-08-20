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

    // Domain change detection in a parallel interpreter run (OpenSeesMP).
    // The Domain stamp is LOCAL to a process, but domainChanged() is COLLECTIVE:
    // the DOF_Numberer numbers the equations of the whole model. So either all
    // the processes call domainChanged() or none of them does; a process that
    // skipped it would also keep equation numbers that a change made in another
    // process has shifted. This returns the logical OR of the flag over all the
    // processes, and the flag itself when the analysis is not run in parallel.
    bool anyDomainChange(bool domainChangedHere);

    // Step-outcome agreement in a parallel interpreter run (OpenSeesMP).
    //
    // Every phase of a step - analysisStep(), newStep(), solveCurrentStep(),
    // commit() - can fail on ONE process: the algorithm does not converge on that
    // subdomain, a material or an element fails its own integration, an
    // integration-error control rejects the step. The analyze() loop then does
    // `return -3` on that process alone. The others know nothing and walk into the
    // next step, whose solve() and getB() are MPI collectives with a participant
    // missing: the job hangs, and if it does not hang it is because the
    // collectives paired up wrongly, which is worse. Verified with
    // scratchpad/asym_exit.tcl - a rank leaving the loop early hangs the rest, and
    // the rank-0 star and the MPI_Allreduce forms of getB() hang identically, so
    // this is not something a cheaper reduction fixes.
    //
    // Returns the WORST result over all the processes (MPI_MINLOC: failures are
    // negative, success is 0), so every process reaches the same verdict and either
    // all of them abandon the step or none does.
    //
    // MUST be called unconditionally at every phase boundary, not only where the
    // local result is negative. A process that called it only on failure would
    // enter this Allreduce while another was inside the solve's Allreduce, and two
    // mismatched collectives on the same communicator is undefined behaviour, not a
    // deadlock one can debug.
    //
    // `phase` names the phase in the diagnostic and is optional; pass a literal.
    int worstStepResult(int resultHere, const char *phase = 0);

    // The two halves of a failure report in OpenSeesMP.
    //
    // A process that failed knows WHY - its element, its material, its
    // convergence history - and says so itself; those messages must never be
    // gated, they are the only ones that carry local detail. What no process
    // knows on its own is WHO failed, and after worstStepResult() every process
    // holds the same verdict, so the N of them print the same line about a
    // failure that happened on one. Hence: local detail from the process that
    // has it, the global verdict once, from process 0.
    //
    // reportHere() is true on process 0 when np > 1, and true always otherwise
    // (serial build, serial run) - so gating a message with it changes nothing
    // outside OpenSeesMP. whoFailed() is the attribution of the last
    // worstStepResult(): worst code, its owner, how many processes failed. It is
    // an empty string when nothing failed and in every serial run, so appending
    // it to a message costs nothing there.
    bool reportHere(void) const;
    const char *whoFailed(void) const;

  private:
    Domain *theDomain;

    // written by worstStepResult(), read by whoFailed(); bounded by construction,
    // the format below truncates the phase name
    char failInfo[160];
};

#endif

