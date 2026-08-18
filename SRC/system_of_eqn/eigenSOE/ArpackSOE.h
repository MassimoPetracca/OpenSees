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
                                                                        
// $Revision: 1.2 $
// $Date: 2009-05-14 22:46:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/ArpackSOE.h,v $

// Written: fmk
// Created: 05/09
//
// Description: This file contains the class definition for ArpackSOE


#ifndef ArpackSOE_h
#define ArpackSOE_h

#include "eigenSOE/EigenSOE.h"
#include <Vector.h>

class AnalysisModel;
class ArpackSolver;
class LinearSOE;
class Channel;

class ArpackSOE : public EigenSOE
{
  public:
    ArpackSOE(double shift = 0.0);

    ~ArpackSOE();

    int setLinks(AnalysisModel &theModel);   
    int setLinearSOE(LinearSOE &theSOE);    

    int getNumEqn(void) const;
    int setSize(Graph &theGraph);
    
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addM(const Matrix &, const ID &, double fact = 1.0);    
   
    void zeroA(void);
    void zeroM(void);

    double getShift(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // wiring for the multi-interpreter (OpenSeesMP) world: set the rank and
    // the channels so that setSize/myMv/checkSameInt take their parallel
    // branches. Mirrors ParallelNumberer/MumpsParallelSOE::setProcessID/
    // setChannels. Without these, processID stays -1 (serial behavior).
    int setProcessID(int processTag);
    int setChannels(int numChannels, Channel **theChannels);

    friend class ArpackSolver;

	int checkSameInt(int);

  protected:
    
  private:
    double *M;
    int Msize;
    bool mDiagonal;

    // C3: the assembled sparse mass matrix.
    //
    // When M is not diagonal, ArpackSolver::myMv used to recompute the whole
    // element loop - elePtr->getM_Force(x) for every FE_Element - once per Lanczos
    // step. Measured on the N=24 cube (45000 DOF, 20 modes, np=2) that is 4.39 s of
    // a 10.6 s eigen loop, 41%: more than the backsolves (3.09 s) and more than the
    // factorisation (2.72 s). It is the same product every time, only the vector
    // changes, so assembling M once turns 279 element loops into 279 sparse
    // matvecs.
    //
    // Compressed column, built from the same Graph that sizes the SOE, so the
    // pattern is guaranteed to hold every (i,j) addM writes - addM already calls
    // addA with the same ID and that succeeds. FULL pattern, both triangles: the
    // mass matrix is symmetric and half of it would do, but a plain matvec has no
    // scatter-writes and this is the first version. Columns of equations this rank
    // does not own stay empty, which keeps M the rank-local PARTIAL that the
    // unconditional reduction in myMv sums - the contract is preserved by
    // construction, see the comment there.
    //
    // Costs Mnnz*(8+4) + (Msize+1)*4 bytes a rank, about +50% on this SOE's matrix
    // memory. Gated on OPS_EIGEN_MSPARSE while the trade is being measured.
    int    *Mcolstart;
    int    *Mrow;
    double *Mval;
    int     Mnnz;
    bool    Msparse;

    int  buildSparseM(Graph &theGraph, int size);
    void freeSparseM(void);
    // slot of (row, col) in the compressed pattern, -1 if the pattern has no such
    // entry (which addM reports rather than silently dropping the term)
    int  sparseMposition(int row, int col) const;
    double shift;
    AnalysisModel *theModel;
    LinearSOE *theSOE;

    int processID;
    int numChannels;
    Channel **theChannels;
    ID **localCol;
    ID *sizeLocal;
};


#endif



