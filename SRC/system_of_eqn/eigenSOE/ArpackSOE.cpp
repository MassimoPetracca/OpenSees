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
                                                                        
// $Revision: 1.3 $
// $Date: 2009-05-14 23:25:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/ArpackSOE.cpp,v $

// Written: fmk
// Created: 05/09
//
// Description: This file contains the class definition for ArpackSOE

#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif

#include <stdlib.h>
#include <ArpackSOE.h>
#include <ArpackSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <f2c.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>



ArpackSOE::ArpackSOE(double s)
:EigenSOE(EigenSOE_TAGS_ArpackSOE),
 M(0), Msize(0), mDiagonal(false),
 Mcolstart(0), Mrow(0), Mval(0), Mnnz(0), Msparse(false),
 shift(s), theModel(0), theSOE(0),
 processID(-1), numChannels(0), theChannels(0), localCol(0), sizeLocal(0)
{
  ArpackSolver *theSolvr = new ArpackSolver();
  this->setSolver(*theSolvr);
  theSolvr->setEigenSOE(*this);
}


int
ArpackSOE::getNumEqn(void) const
{
  if (theSOE != 0)
    return theSOE->getNumEqn();
  else 
    return 0;
}
    
ArpackSOE::~ArpackSOE()
{
  if (M != 0) delete [] M;
  this->freeSparseM();
  if (theChannels != 0) delete [] theChannels;
}

// Assemble M once instead of redoing the element loop per Lanczos step. This is a
// FEATURE switch, not a diagnostic: it buys latency and pays memory, so which way it
// should default is a modelling decision. Measured on the N=24 cube (45000 DOF, 20
// modes): the mass product goes 10x faster and the whole modal analysis 1.3-2.3x,
// against Mnnz*(8+4) + (Msize+1)*4 extra bytes a rank - about +50% on this SOE's
// matrix memory, ~20 MB at 45000 DOF but ~3 GB at n=1e7.
static bool opsEigenMsparseOn(void) {
  static const bool on = (getenv("OPS_EIGEN_MSPARSE") != 0);
  return on;
}

void
ArpackSOE::freeSparseM(void)
{
  if (Mcolstart != 0) { delete [] Mcolstart; Mcolstart = 0; }
  if (Mrow != 0)      { delete [] Mrow;      Mrow = 0; }
  if (Mval != 0)      { delete [] Mval;      Mval = 0; }
  Mnnz = 0;
  Msparse = false;
}

int
ArpackSOE::sparseMposition(int row, int col) const
{
  // rows are stored ascending within a column, so the scan can stop early - the
  // same shape as MumpsSOE::addA, whose cost this mirrors
  for (int k = Mcolstart[col]; k < Mcolstart[col+1]; k++) {
    if (Mrow[k] == row)
      return k;
    if (Mrow[k] > row)
      return -1;
  }
  return -1;
}

int
ArpackSOE::buildSparseM(Graph &theGraph, int size)
{
  this->freeSparseM();

  if (size <= 0)
    return 0;

  // Pass 1: one slot for the diagonal plus one per adjacency entry. Columns of
  // equations this rank does not own have no vertex in the local graph and stay
  // empty - which is what keeps M a per-rank PARTIAL.
  Mcolstart = new int[size+1];
  if (Mcolstart == 0) {
    opserr << "ArpackSOE::buildSparseM - out of memory for the column pointers\n";
    return -1;
  }

  Mcolstart[0] = 0;
  for (int a = 0; a < size; a++) {
    int len = 0;
    Vertex *theVertex = theGraph.getVertexPtr(a);
    if (theVertex != 0)
      len = 1 + theVertex->getAdjacency().Size();
    Mcolstart[a+1] = Mcolstart[a] + len;
  }
  Mnnz = Mcolstart[size];

  Mrow = new int[(Mnnz > 0) ? Mnnz : 1];
  Mval = new double[(Mnnz > 0) ? Mnnz : 1];
  if (Mrow == 0 || Mval == 0) {
    opserr << "ArpackSOE::buildSparseM - out of memory for " << Mnnz
	   << " nonzeros\n";
    this->freeSparseM();
    return -1;
  }

  // Pass 2: fill the rows, kept ascending by insertion so sparseMposition() can
  // stop early. Columns are short (mean 36 on the N=24 cube), so insertion is the
  // right sort here.
  for (int a = 0; a < size; a++) {
    Vertex *theVertex = theGraph.getVertexPtr(a);
    if (theVertex == 0)
      continue;
    int lo = Mcolstart[a];
    int at = lo;
    Mrow[at++] = theVertex->getTag();       // diagonal first, then insert the rest
    const ID &adj = theVertex->getAdjacency();
    int adjSize = adj.Size();
    for (int i = 0; i < adjSize; i++) {
      int row = adj(i);
      int j = at - 1;
      while (j >= lo && Mrow[j] > row) {
	Mrow[j+1] = Mrow[j];
	j--;
      }
      Mrow[j+1] = row;
      at++;
    }
  }

  for (int k = 0; k < Mnnz; k++)
    Mval[k] = 0.0;

  Msparse = true;
  return 0;
}

int
ArpackSOE::setProcessID(int processTag)
{
  processID = processTag;
  return 0;
}

int
ArpackSOE::setChannels(int nChannels, Channel **theC)
{
  numChannels = nChannels;

  if (theChannels != 0)
    delete [] theChannels;

  theChannels = new Channel *[numChannels];
  for (int i=0; i<numChannels; i++)
    theChannels[i] = theC[i];

  return 0;
}

int 
ArpackSOE::setSize(Graph &theGraph)
{
  if (theSOE == 0)
    return -1;

  int result = 0;
  int size = 0;

  if (processID == -1) {

    size = theGraph.getNumVertex();    

  } else {

    // fist itearte through the vertices of the graph to get n
    int maxVertexTag = -1;
    Vertex *theVertex;
    
    VertexIter &theVertices = theGraph.getVertices();
    while ((theVertex = theVertices()) != 0) {
      int vertexTag = theVertex->getTag();
      if (vertexTag > maxVertexTag)
	maxVertexTag = vertexTag;
    }

    if (processID != 0) {
      
      //
      // if subprocess, send local max vertexTag (n)
      // recv ax n from P0
      //

      static ID data(1);
      
      data(0) = maxVertexTag;
      Channel *theChannel = theChannels[0];
      theChannel->sendID(0, 0, data);    // send local max tag
      theChannel->recvID(0, 0, data);    // recv global max tag 

      size = data(0);    
    } else {

      //
      // from each distributed soe recv it's max n and compare; return max n to all
      //
      
      static ID data(1);
      FEM_ObjectBroker theBroker;

      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->recvID(0, 0, data);
	if (data(0) > maxVertexTag)
	  maxVertexTag = data(0);
      }
      
      data(0) = maxVertexTag;

      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->sendID(0, 0, data);
      }
      size = maxVertexTag;
    }

    size+=1; // vertices numbered 0 through n-1
  }
  
  // invoke setSize on the linearSOE
  /* SETSIZE ALREADY CALLED IF USING ARPACK AS THIS SOE DOING DOUBLE DUTY
  result = theSOE->setSize(theGraph);
  if (result < 0) {
    opserr << "WARNING ArpackSOE::ArpackSOE : - LinearSOE - failed in setSize\n";
    return -1;
  }
  */

  if (size != Msize && size > 0) {

    if (M != 0) 
      delete [] M;
    
    M = new double[size];
    
    if (M == 0) {
      opserr << "WARNING ArpackSOE::ArpackSOE : - out of memory creating memory for M\n";
      Msize = 0;
    } else
      Msize = size;
  }

  // C3: the sparse pattern has to be rebuilt whenever the graph changes, not only
  // when the size changes - a domain change can keep the equation count and move
  // the connectivity. Cheap relative to what it saves, and it is per setSize, not
  // per solve. Failure is not fatal: Msparse stays false and myMv keeps the element
  // loop.
  if (opsEigenMsparseOn()) {
    if (this->buildSparseM(theGraph, size) < 0)
      opserr << "WARNING ArpackSOE::setSize - sparse M unavailable, the mass "
	     << "product falls back to the element loop\n";
  } else
    this->freeSparseM();

  //
  // invoke setSize() on the Solver
  //

  EigenSolver *theSolvr = this->getSolver();

  if (theSolvr == 0) {
    opserr << "ArpackSOE::setSize(Graph &theGraph) - no EigenSolver set\n";             
    return -1;
  }
  int solverOK = theSolvr->setSize();

  if (solverOK < 0) {
    opserr << "WARNING:ArpackSOE::setSize() -  solver failed setSize()\n";
    return solverOK;
  } 
  
  return result;    
}

int 
ArpackSOE::addA(const Matrix &m, const ID &id, double fact)
{
  if (theSOE == 0) {
    opserr << "ArpackSOE::addA() - no SOE set\n";
    return -1;
  }

  // check for a quick return 
  if (fact == 0.0)  return 0;

  return theSOE->addA(m, id, fact);
}


void 
ArpackSOE::zeroA(void)
{
  if (theSOE == 0) {
    opserr << "ArpackSOE::zeroA() - no SOE set\n";
    return;
  }
  return theSOE->zeroA();
}

int 
ArpackSOE::addM(const Matrix &m, const ID &id, double fact)
{
  if (theSOE == 0) {
    opserr << "ArpackSOE::addM() - no SOE set\n";
    return -1;
  }

  int res = this->addA(m, id, -shift);

  if (res < 0)
    return res;

  int idSize = id.Size();

  // C3: scatter into the assembled sparse M. Unconditional - it must run for every
  // contribution, including those after mDiagonal has been cleared, which is where
  // the diagonal bookkeeping below stops. Only the entries whose BOTH indices are
  // in range are taken, exactly like the diagonal path: a negative equation number
  // is a constrained dof and carries no mass term.
  if (Msparse) {
    for (int i=0; i<idSize; i++) {
      int locI = id(i);
      if (locI < 0 || locI >= Msize)
	continue;
      for (int j=0; j<idSize; j++) {
	int locJ = id(j);
	if (locJ < 0 || locJ >= Msize)
	  continue;
	if (m(i,j) == 0.0)
	  continue;
	int k = this->sparseMposition(locI, locJ);
	if (k >= 0)
	  Mval[k] += m(i,j);
	else {
	  // The graph that built this pattern is the same one addA works from, so
	  // a miss means the two have gone out of step. Report it - silently
	  // dropping a mass term would move the eigenvalues and look plausible.
	  opserr << "ArpackSOE::addM - no slot for (" << locI << "," << locJ
		 << ") in the sparse mass pattern; disabling it and falling back "
		 << "to the element loop\n";
	  this->freeSparseM();
	  break;
	}
      }
      if (Msparse == false)
	break;
    }
  }

  if (mDiagonal == false)
    return  res;
  for (int i=0; i<idSize; i++) {
    int locI = id(i);
    if (locI >= 0 && locI < Msize) {
      for (int j=0; j<idSize; j++) {
	int locJ = id(j);
	if (locJ >= 0 && locJ < Msize) {
	  if (locI == locJ) {
	    // m(i,j) not m(i,i): duplicate IDs (Transformation) must sum the block.
	    M[locI] += m(i,j);
	  } else {
	    if (m(i,j) != 0.0) {
	      mDiagonal = false;
	      return res;
	    }
	  }
	}
      }
    }
  }

  return 0;
}   
 
void 
ArpackSOE::zeroM(void)
{
  if (theSOE == 0) {
    opserr << "ArpackSOE::zeroM() - no SOE set\n";
    return;
  }

  mDiagonal = true;

  for (int i=0; i<Msize; i++)
    M[i] = 0;

  if (Msparse)
    for (int k=0; k<Mnnz; k++)
      Mval[k] = 0.0;
}


double 
ArpackSOE::getShift(void)
{
    return shift;
}


int 
ArpackSOE::sendSelf(int commitTag, Channel &theChannel)
{
  int sendID =0;
  
  if (processID == -1)
    processID = 0;
  
  // if P0 check if already sent. If already sent use old processID; if not allocate a new process 
  // id for remote part of object, enlarge channel * to hold a channel * for this remote object.
  // if not P0, send current processID

  if (processID == 0) {
    
    // check if already using this object
    bool found = false;
    for (int i=0; i<numChannels; i++)
      if (theChannels[i] == &theChannel) {
	sendID = i+1;
	found = true;
      }
    
    // if new object, enlarge Channel pointers to hold new channel * & allocate new ID
    if (found == false) {
      int nextNumChannels = numChannels + 1;
      Channel **nextChannels = new Channel *[nextNumChannels];
      if (nextNumChannels == 0) {
	opserr << "ArpackSOE::sendSelf() - failed to allocate channel array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	nextChannels[i] = theChannels[i];
      nextChannels[numChannels] = &theChannel;
      
      numChannels = nextNumChannels;
      
      if (theChannels != 0)
	delete [] theChannels;
      
      theChannels = nextChannels;
      
      if (localCol != 0)
	delete [] localCol;
      localCol = new ID *[numChannels];
      if (localCol == 0) {
	opserr << "ArpackSOE::sendSelf() - failed to allocate id array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	localCol[i] = 0;    
      
      if (sizeLocal != 0)
	delete sizeLocal;
      
      sizeLocal = new ID(numChannels);
      
      // allocate new processID for remote object
      sendID = numChannels;
    }
  } else 
    sendID = processID;
  
  // send remotes processID
  ID idData(1);
  idData(0) = sendID;

  int res = theChannel.sendID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING ArpackSOE::sendSelf() - failed to send data\n";
    return -1;
  }

  return 0;  
}

    
int 
ArpackSOE::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  ID idData(1);
  int res = theChannel.recvID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING ArpackSOE::recvSelf() - failed to send data\n";
    return -1;
  }	      
  processID = idData(0);

  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  localCol = new ID *[numChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;

  if (sizeLocal != 0)
    delete sizeLocal;

  sizeLocal = new ID(numChannels);

  return 0;
}

int 
ArpackSOE::setLinks(AnalysisModel &theAnalysisModel)
{
  theModel = &theAnalysisModel;
  return 0;
}

int 
ArpackSOE::setLinearSOE(LinearSOE &theLinearSOE)
{
  theSOE = &theLinearSOE;
  return 0;
}

int
ArpackSOE::checkSameInt(int value)
{
	if (processID == -1)
		return 1;

#ifdef _PARALLEL_INTERPRETERS

	// The last rank-0 star on the eigen path, and the one that runs most often:
	// ArpackSolver calls this once per Lanczos iteration to check that every rank
	// got the same `ido` out of its own dsaupd (368 times on the N=24 cube, against
	// 279 mass products and 88 backsolves). Two ints and one collective replace
	// O(P) messages through a serial owner, and every rank ends with the same
	// verdict by construction rather than because rank 0 posted it back.
	//
	// max and min in one reduction: d[0] carries the value, d[1] its negation, so
	// MPI_MAX gives max in d[0] and -min in d[1]. They agree iff every rank passed
	// the same value.
	int d[2];
	d[0] = value;
	d[1] = -value;
	if (MPI_Allreduce(MPI_IN_PLACE, d, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD)
	    != MPI_SUCCESS) {
		opserr << "ArpackSOE::checkSameInt() - MPI_Allreduce failed\n";
		return 0;
	}
	return (d[0] == -d[1]) ? 1 : 0;

#else

	static ID idData(1);
    if (processID != 0) {

		Channel *theChannel = theChannels[0];
	    idData(0) = value;
		theChannel->sendID(0, 0, idData);
		theChannel->recvID(0, 0, idData);
		if (idData(0) == 1)
			return 1;
		else
			return 0;
	}

	else {
        int ok = 1;
		// receive B 
		for (int j=0; j<numChannels; j++) {
		// get X & add
			Channel *theChannel = theChannels[j];
			theChannel->recvID(0, 0, idData);
			if (idData(0) != value)
				ok = 0;
		}

		// send results back
		idData(0) = ok;
		for (int j=0; j<numChannels; j++) {
			Channel *theChannel = theChannels[j];
			theChannel->sendID(0, 0, idData);
		}
		return ok;
    }

#endif
}
