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
                                                                        
// $Revision: 1.5 $
// $Date: 2009-05-11 20:56:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.h,v $
                                                                        
#ifndef MumpsParallelSOE_h
#define MumpsParallelSOE_h

// Written: fmk 
// Description: This file contains the class definition for MumpsParallelSOE
// MumpsParallelSOE is a subclass of LinearSOE. It uses the sparse column
// storage scheme. The matrix A is kept distributed, X and B kept on all processors.
//
// matrix types (matType): 0 Unsymmetrc
//                         1 Symmetrix positive definite
//                         2 General Symmetric

// What: "@(#) MumpsParallelSOE.h, revA"


#include <mpi.h>

#include <MumpsSOE.h>
#include <Vector.h>

class MumpsParallelSolver;

class MumpsParallelSOE : public MumpsSOE
{
  public:
    MumpsParallelSOE(MumpsParallelSolver &theSolver, int matType=0);
    MumpsParallelSOE();
    
    ~MumpsParallelSOE();

    // these methods need to be rewritten
    int setSize(Graph &theGraph);

    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);            
    const Vector &getB(void);
    void zeroB(void);
    int solve(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    friend class MumpsParallelSolver;        

    int setProcessID(int processTag);
    int setChannels(int numChannels, Channel **theChannels);

  protected:
    // Does B still hold the sum of every rank's myB?
    //
    // getB() and solve() both merge myB into B, and in a Newton iteration they do
    // it one after the other on data only one of them changed: solve() merges the
    // residual the previous formUnbalance() built, then formUnbalance() rebuilds
    // myB, then the convergence test merges it again through getB(). Whichever
    // runs second is re-summing what the first already summed. This flag lets the
    // second one skip: raised by every writer of myB (addB, setB, zeroB, setSize)
    // and cleared by whoever merges.
    //
    // The decision to merge is taken COLLECTIVELY, not from this flag alone. The
    // merge is an MPI collective, so a rank that skipped it while another entered
    // it would hang the job - the same asymmetry that used to deadlock the staged
    // analysis. An MPI_Allreduce of one int costs microseconds against the
    // size*8 bytes the merge moves, so agreeing is far cheaper than the work it
    // saves.
    //
    // Invariant this rests on: B is read through getB() and nowhere else.
    // MumpsSOE::normRHS() does read B directly, and it has no caller anywhere in
    // SRC (checked); it is deliberately NOT made to merge, because that would put
    // a collective on a path nothing drives in lockstep. If it ever gains a
    // caller, it must call getB() first.
    bool myBdirty;

    // Does B hold the global sum on EVERY rank, or only on rank 0?
    //
    // The two readers of B want different things. MUMPS reads it only on the host
    // (MumpsParallelSolver copies B into X and hands id.rhs over on rank 0 alone),
    // so solve() only owes the sum to rank 0 - an MPI_Reduce, half the volume of an
    // Allreduce and no return leg. getB() hands the vector to the convergence
    // tests on every rank and does owe the full Allreduce. Recording how far the
    // last merge got is what lets solve() take the cheap one without leaving a
    // later getB() reading a stale B on the other ranks.
    //
    // In a modal analysis this is the whole RHS exchange: ArpackSolver installs the
    // global vector with setB() on rank 0 and calls zeroB() on the others, so the
    // sum IS rank 0's own myB and no communication is needed at all. That happens
    // once per backsolve, hundreds of times per analysis.
    bool Bglobal;

    // How far the pending merge of myB into B has to go, for a caller that needs B
    // on every rank (needGlobal true, i.e. getB) or on the host alone (needGlobal
    // false, i.e. solve). Collective: every rank must come out with the same
    // answer, or one of them is left alone in a reduction. One Allreduce of four
    // ints decides it, and that includes agreeing on Bglobal itself rather than
    // trusting the local copy of the flag.
    //
    //   MERGE_NONE     B is already as far along as this caller needs
    //   MERGE_LOCAL    rank 0's own myB is the whole sum: copy it, no communication
    //   MERGE_REDUCE   sum into rank 0 only (MPI_Reduce)
    //   MERGE_BCAST    rank 0 already holds the sum, hand it to the others
    //   MERGE_ALL      sum into every rank (MPI_Allreduce)
    enum { MERGE_NONE = 0, MERGE_LOCAL, MERGE_REDUCE, MERGE_BCAST, MERGE_ALL };
    int planMergeOfB(bool needGlobal);



    
  private:
    int processID;
    int numChannels;
    Channel **theChannels;
    ID **localCol;

    // receive buffer of the rank-0 star in getB(); allocated only where that star
    // is compiled in, i.e. NOT under _PARALLEL_INTERPRETERS. The old sizeWork
    // companion was initialised to 0 and never assigned again - it is gone.
    double *workArea;
    double *myB;
    Vector *myVectB;

    // set by setSize(): the reduced global equation count is zero, i.e. every
    // DOF of the model is prescribed. Distinct from "size is still 0 because
    // setSize() was never called on this rank", which must NOT skip the solve
    // (the other ranks would block in the collectives).
    bool zeroEqnSystem;
};


#endif

