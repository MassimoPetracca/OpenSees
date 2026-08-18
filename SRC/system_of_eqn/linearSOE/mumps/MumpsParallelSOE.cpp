/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
* University of California, Berkeley, is strictly prohibited.  See   **
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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp,v $
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the implementation for MumpsParallelSOE

#include <stdlib.h>
#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif

#include <MumpsParallelSOE.h>
#include <MumpsParallelSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


// reported once per process: a system with no unknowns is a legitimate state,
// not an error, but it should be visible in the log
static bool opsMumpsZeroEqnNoted = false;

MumpsParallelSOE::MumpsParallelSOE(MumpsParallelSolver &theSolvr, int matType)
  :MumpsSOE(theSolvr, LinSOE_TAGS_MumpsParallelSOE, matType),
   processID(0), numChannels(0), theChannels(0), localCol(0), workArea(0),
   myB(0), myVectB(0), zeroEqnSystem(false), myBdirty(true), Bglobal(false)
{
    theSolvr.setLinearSOE(*this);
}


MumpsParallelSOE::MumpsParallelSOE()
  :MumpsSOE(LinSOE_TAGS_MumpsParallelSOE),
   processID(0), numChannels(0), theChannels(0), localCol(0), workArea(0),
   myB(0), myVectB(0), zeroEqnSystem(false), myBdirty(true), Bglobal(false)
{

}


MumpsParallelSOE::~MumpsParallelSOE()
{

  if (theChannels != 0)
    delete [] theChannels;

  if (localCol != 0)
    for (int i=0; i<numChannels; i++)
      if (localCol[i] != 0)
	delete localCol[i];
  delete [] localCol;

  if (myB != 0)
    delete [] myB;

  if (myVectB != 0)
    delete myVectB;
}

int 
MumpsParallelSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int oldSize = size;
  int maxNumSubVertex = 0;
  
  // fist itearte through the vertices of the graph to get nnzLoc and n
  int maxVertexTag = -1;
  Vertex *theVertex;
  int newNNZ = 0;
  size = theGraph.getNumVertex();
  int mySize = size;
  //opserr << "MumpsParallelSOE: size : " << size << endln;

  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    int vertexTag = theVertex->getTag();
    if (vertexTag > maxVertexTag)
      maxVertexTag = vertexTag;
    const ID &theAdjacency = theVertex->getAdjacency();
    newNNZ += theAdjacency.Size() +1; // the +1 is for the diag entry
  }

  if (matType !=  0) {

    // symmetric - allows us to reduce nnz by almost half
    newNNZ -= size;
    newNNZ /= 2;
    newNNZ += size;
  }

  nnz = newNNZ;

  if (processID != 0) {

    //
    // if subprocess, send local max vertexTag (n)
    // recv ax n from P0
    //
    static ID data(1);

    data(0) = maxVertexTag;
    Channel *theChannel = theChannels[0];
    theChannel->sendID(0, 0, data);
    theChannel->recvID(0, 0, data);
    
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

  if (nnz > Asize) { // we have to get more space for A and rowA and colA

    if (A != 0) delete [] A;
    if (rowA != 0) delete [] rowA;
    if (colA != 0) delete [] colA;

    A = new double[nnz];    
    rowA = new int[nnz];
    colA = new int[nnz];
      
    for (int i=0; i<nnz; i++) {
      A[i]=0;
      rowA[i]=0;
      colA[i]=0;
    }

    if (rowA == 0 || A == 0 || colA == 0) {
      opserr << "WARNING SparseGenColLinSOE::SparseGenColLinSOE :";
      opserr << " ran out of memory for A and rowA with nnz = ";
      opserr << nnz << " \n";
      size = 0; Asize = 0; nnz = 0;
	result =  -1;
    } 
    Asize = nnz;
  }

  if (size > Bsize) { // we have to get space for the vectors

    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (myB != 0) delete [] myB;
    if (workArea != 0) delete [] workArea;
    if (colStartA != 0)  delete [] colStartA;

    // create the new
    B = new double[size];
    X = new double[size];
    myB = new double[size];
    colStartA = new int[size+1];

    // workArea is the receive buffer of the rank-0 star in getB(), and that star
    // now exists only under _PARALLEL_PROCESSING: in an OpenSeesMP build getB()
    // reduces with MPI_Allreduce and nothing ever reads workArea. Allocating it
    // there costs size doubles per rank - the same as B, X and myB, so a quarter
    // of this SOE's vector memory - for nothing. At N=64 that is 6.5 MB a rank,
    // and it is per rank, not per job.
#ifndef _PARALLEL_INTERPRETERS
    workArea = new double[size];
    if (workArea == 0) {
      opserr << "WARNING MumpsParallelSOE::setSize - ran out of memory for the"
	     << " work area (size " << size << ")\n";
      size = 0; Bsize = 0;
      return -1;
    }
#endif

    if (B == 0 || X == 0 || colStartA == 0 || myB == 0) {
      opserr << "WARNING MumpsSOE::MumpsSOE :";
      opserr << " ran out of memory for vectors (size) (";
      opserr << size << ") \n";
      size = 0; Bsize = 0;
      result =  -1;
    }
    else
      Bsize = size;

  }
  
  // zero the vectors
  for (int j=0; j<size; j++) {
    B[j] = 0;
    X[j] = 0;
    myB[j] = 0;
  }

  // B and myB are both zero here, so B does hold their sum - but the arrays may
  // also have just been reallocated, so claim nothing and let the next consumer
  // merge once
  myBdirty = true;
  Bglobal = false;

  // create new Vectors objects. The (vectX == 0) terms matter when size is 0:
  // size == oldSize == 0 would otherwise leave them unallocated, and
  // getX()/getB() are FATAL on a null pointer.
  if (size != oldSize || vectX == 0 || vectB == 0 || myVectB == 0) {
    if (vectX != 0) delete vectX;
    if (vectB != 0) delete vectB;
    if (myVectB != 0) delete myVectB;

    vectX = new Vector(X,size);
    vectB = new Vector(B,size);
    myVectB = new Vector(myB, size);
  }

  //
  // a system with no unknowns
  //
  // Every DOF of the model is prescribed, so the constraint handler has
  // eliminated all of them: a staged analysis step in which nothing is active,
  // or a displacement-controlled single-element test. The system A x = b has no
  // unknowns and its unique solution is the empty vector - there is nothing to
  // factor and nothing to solve. The prescribed values themselves are enforced
  // by the constraint handler (applyLoad -> enforceSPs), not by this solve, so
  // the step still produces the correct displacements, element states and
  // reactions.
  //
  // MUMPS must not be entered at all here: id.n = 0 returns INFO(1) = -16
  // ("N out of range"). size is the reduced global equation count and is
  // therefore identical on every rank, so this branch is taken collectively and
  // the ranks stay in step.
  //
  zeroEqnSystem = (size == 0);
  if (zeroEqnSystem) {
    factored = false;
    if (opsMumpsZeroEqnNoted == false) {
      opsMumpsZeroEqnNoted = true;
      opserr << "MumpsParallelSOE::setSize - the model has no free equations: "
	     << "all DOFs are prescribed by the constraint handler. The linear "
	     << "solve is skipped; prescribed values are still enforced.\n";
    }
    return result;
  }

  // fill in colStartA and rowA
  if (size != 0) {
    colStartA[0] = 0;
    int startLoc = 0;
    int lastLoc = 0;
    for (int a=0; a<size; a++) {
      
      theVertex = theGraph.getVertexPtr(a);
      if (theVertex != 0) {
	
	int vertexTag = theVertex->getTag();
	rowA[lastLoc++] = vertexTag; // place diag in first
	const ID &theAdjacency = theVertex->getAdjacency();
	int idSize = theAdjacency.Size();
	
	// now we have to place the entries in the ID into order in rowA
	
	if (matType != 0) {
	  
	  // symmetric
	  for (int i=0; i<idSize; i++) {
	    int row = theAdjacency(i);
	    if (row > vertexTag) {
	      bool foundPlace = false;
	      // find a place in rowA for current col
	      for (int j=startLoc; j<lastLoc; j++)
		if (rowA[j] > row) { 
		  // move the entries already there one further on
		  // and place col in current location
		  for (int k=lastLoc; k>j; k--)
		    rowA[k] = rowA[k-1];
		  rowA[j] = row;
		  foundPlace = true;
		  j = lastLoc;
		}
	      
	      if (foundPlace == false) // put in at the end
		rowA[lastLoc] = row;
	      lastLoc++;
	    }
	  }
	  
	} else {

	  // unsymmetric	  
	  for (int i=0; i<idSize; i++) {
	    int row = theAdjacency(i);
	    bool foundPlace = false;
	    // find a place in rowA for current col
	    for (int j=startLoc; j<lastLoc; j++)
	      if (rowA[j] > row) { 
		// move the entries already there one further on
		// and place col in current location
		for (int k=lastLoc; k>j; k--)
		  rowA[k] = rowA[k-1];
		rowA[j] = row;
		foundPlace = true;
		j = lastLoc;
	      }
	    if (foundPlace == false) // put in at the end
	      rowA[lastLoc] = row;
	    
	    lastLoc++;
	  }
	}
      }
      colStartA[a+1] = lastLoc;
      startLoc = lastLoc;
    }
  }

  // fill in colA
  int count = 0;
  for (int i=0; i<size; i++) {
    for (int k=colStartA[i]; k<colStartA[i+1]; k++)
      colA[count++] = i;
  }

  // this rank builds its own rowA above, so it owes addA() the same ascending
  // invariant the base class checks - the columns of the dofs this rank does not
  // own stay empty, which the merge in addA() handles as an immediate miss
  if (this->verifyRowAOrder() < 0)
    return -1;

  LinearSOESolver *theSolvr = this->getSolver();

  int solverOK = theSolvr->setSize();
  if (solverOK < 0) {
    opserr << "WARNING:MumpsParallelSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    

  return result;    
}


int 
MumpsParallelSOE::solve(void)
{
  int resSolver = 0;

  // no unknowns: the empty vector is the solution, see setSize(). Keyed on
  // zeroEqnSystem, not on size: this rank must have been through setSize(), so
  // that every rank agrees and none is left waiting in a collective.
  if (zeroEqnSystem)
    return 0;

  //
  // if subprocess send B, solve and recv back X and B
  //


#ifdef _PARALLEL_INTERPRETERS

  // OpenSeesMP: solve() is entered by every rank in lockstep and the MUMPS
  // job=3 below is collective on MPI_COMM_WORLD, so the star exchange over
  // Channels can be replaced by two tree collectives with the same final
  // state on every rank: B = sum of the locally assembled myB (the sum is
  // exactly what the old rank-0 gather computed and shipped back), X = the
  // solution broadcast from the host.

  // ... and only when myB has moved since the last merge. In a Newton iteration
  // the convergence test has usually just merged this same myB through getB(),
  // and re-summing it would be a second full exchange of identical data.
  //
  // And only as far as the host. MumpsParallelSolver reads B on rank 0 alone -
  // "if (rank == 0) { X[i] = B[i]; id.rhs = X; }" - so what solve() owes is the sum
  // ON RANK 0, not on every rank. MPI_Reduce instead of MPI_Allreduce: half the
  // volume and no return leg. Bglobal records that the other ranks were left
  // behind, so a getB() that follows still gets a correct vector (one Bcast, which
  // is the leg the Allreduce was paying for anyway).
  //
  // Better still when one rank holds the whole RHS and it is the host: then there is
  // nothing to communicate at all. That is not a special case, it is the modal
  // analysis - ArpackSolver installs the global vector with setB() on rank 0 and
  // calls zeroB() on the others, once per backsolve.
  int plan = this->planMergeOfB(false);

  if (plan == MERGE_LOCAL) {
    if (processID == 0)
      for (int i = 0; i < size; i++)
	B[i] = myB[i];
    myBdirty = false;
    Bglobal = false;
  } else if (plan == MERGE_REDUCE) {
    MPI_Reduce(myB, B, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    myBdirty = false;
    Bglobal = false;
  }


  resSolver = this->LinearSOE::solve();

  if (resSolver == 0) {
    MPI_Bcast(X, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (processID != 0)
      factored = true;
  }

#else

  if (processID != 0) {

    // send B
      Channel *theChannel = theChannels[0];
    theChannel->sendVector(0, 0, *myVectB);
  
    resSolver =  this->LinearSOE::solve();
  
    if (resSolver == 0) {
      // receive X,B and result
      theChannel->recvVector(0, 0, *vectX);
      theChannel->recvVector(0, 0, *vectB);
      factored = true;
    }
    }

  //
  // if main process, recv B & A from all, solve and send back X, B & result
  //

  else {

    // add P0 contribution to B
      *vectB = *myVectB;

    // receive B
    for (int j=0; j<numChannels; j++) {
      // get X & add
      Channel *theChannel = theChannels[j];
      theChannel->recvVector(0, 0, *vectX);
      *vectB += *vectX;
    }
  
    // solve
    resSolver = this->LinearSOE::solve();
  
    // send results back
    if (resSolver == 0) {
      for (int j=0; j<numChannels; j++) {
	Channel *theChannel = theChannels[j];
	theChannel->sendVector(0, 0, *vectX);
	theChannel->sendVector(0, 0, *vectB);
      }
    }
    }

#endif

  return resSolver;
}



int 
MumpsParallelSOE::addB(const Vector &v, const ID &id, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;

  int idSize = id.Size();    
  // check that m and id are of similar size
  if (idSize != v.Size() ) {
    opserr << "SparseGenColLinSOE::addB() ";
    opserr << " - Vector and ID not of similar sizes\n";
    return -1;
  }    

  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	myB[pos] += v(i);
    }
  } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	myB[pos] -= v(i);
    }
  } else {
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
	myB[pos] += v(i) * fact;
    }
  }

  myBdirty = true;
  Bglobal = false;

  return 0;
}


int
MumpsParallelSOE::setB(const Vector &v, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;

  //opserr << "MumpsParallelSOE::setB() - start()\n";
  //opserr << v;

  if (v.Size() != size) {
    opserr << "WARNING MumpsParallelSOE::setB() -";
    opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
    return -1;
  }
    
  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<size; i++) {
      myB[i] = v(i);
    }
  } else if (fact == -1.0) {
    for (int i=0; i<size; i++) {
      myB[i] = -v(i);
    }
  } else {
    for (int i=0; i<size; i++) {
      myB[i] = v(i) * fact;
    }
  }

   //opserr << "MumpsParallelSOE::setB() - end()\n";
  myBdirty = true;
  Bglobal = false;
  return 0;
}

void
MumpsParallelSOE::zeroB(void)
{
  double *Bptr = myB;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;

  myBdirty = true;
  Bglobal = false;
}


int
MumpsParallelSOE::planMergeOfB(bool needGlobal)
{
#ifdef _PARALLEL_INTERPRETERS

  // Four integers, one MPI_SUM, decided collectively - the same reason
  // the #0/#0b lesson: the merge is a collective and a
  // rank that took a different branch would leave the others waiting.
  //
  //   d[0]  has this rank written myB since the last merge?
  //   d[1]  is this rank's myB anything other than identically zero?
  //   d[2]  rank+1 if so, 0 otherwise
  //   d[3]  does this rank believe B is NOT the global sum yet?
  //
  // After the sum: d[1] counts the ranks that actually carry a contribution, and
  // when that count is 1 the value d[2]-1 names the one rank that does. The test is
  // exact - it reads the data, it does not guess - and the scan that feeds it stops
  // at the first nonzero, so it costs O(1) whenever there IS a contribution and
  // O(size) compares only on the all-zero ranks, against the size*8 bytes of
  // network it removes. d[3] is in the same reduction so that Bglobal is agreed
  // rather than assumed: the flag is cleared by addB/setB/zeroB, and nothing
  // guarantees every rank calls one of them in a given iteration.
  int localNonzero = 0;
  if (needGlobal == false) {
    for (int i = 0; i < size; i++)
      if (myB[i] != 0.0) { localNonzero = 1; break; }
  }

  int d[4];
  d[0] = myBdirty ? 1 : 0;
  d[1] = localNonzero;
  d[2] = localNonzero ? (processID + 1) : 0;
  d[3] = Bglobal ? 0 : 1;

  if (MPI_Allreduce(MPI_IN_PLACE, d, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD)
      != MPI_SUCCESS) {
    opserr << "MumpsParallelSOE::planMergeOfB - MPI_Allreduce failed; falling back "
	   << "to the full merge\n";
    return MERGE_ALL;
  }

  const bool anyDirty = (d[0] != 0);
  const bool anyStale = (d[3] != 0);

  if (needGlobal) {
    if (anyDirty)  return MERGE_ALL;
    if (anyStale)  return MERGE_BCAST;   // rank 0 has the sum, the others do not
    return MERGE_NONE;
  }

  if (!anyDirty)
    return MERGE_NONE;                   // B on the host is as merged as it was

  // Only rank 0's copy has to be right, so the cheap case is "rank 0 already holds
  // the whole sum". Two ways that happens: nobody contributes anything (d[1] == 0,
  // the sum is zero and rank 0's zero myB is exactly it) or exactly one rank
  // contributes and it is rank 0 - which is what ArpackSolver does on every
  // backsolve, setB() on the host and zeroB() everywhere else.
  if (d[1] == 0 || (d[1] == 1 && d[2] == 1))
    return MERGE_LOCAL;

  return MERGE_REDUCE;

#else

  return needGlobal ? MERGE_ALL : MERGE_REDUCE;

#endif
}



const Vector &
MumpsParallelSOE::getB(void)
{
  // no unknowns: nothing to merge across the ranks, see setSize()
  if (zeroEqnSystem)
    return *vectB;

#ifdef _PARALLEL_INTERPRETERS

  // Same reasoning as solve(): getB() is entered by every rank in lockstep - it
  // is called by the convergence tests that measure the unbalance, by the
  // accelerators and by the line searches, all driven by the same algorithm on
  // every rank - so the star exchange over the Channels can be replaced by the
  // tree collective that computes the same sum. This is the form
  // ClusterPardisoSOE::getB() already uses.
  //
  // Skipped when nothing has written myB since the last merge: B is then still
  // the sum every rank agreed on. solve() carries the same guard, so between the
  // two of them the RHS crosses the network once per iteration instead of twice.
  //
  // Third case since solve() stopped merging beyond the host: myB is clean but B
  // is the sum on rank 0 only. The callers of getB() are on every rank, so the sum
  // has to travel the last leg - and a Bcast is exactly that leg, half the volume
  // of the Allreduce, because rank 0 already holds the total.
  {
    int plan = this->planMergeOfB(true);

    if (plan == MERGE_ALL) {
      MPI_Allreduce(myB, B, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      myBdirty = false;
      Bglobal = true;
    } else if (plan == MERGE_BCAST) {
      MPI_Bcast(B, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      Bglobal = true;
    }
  }

#else

  if (processID != 0) {
    Channel *theChannel = theChannels[0];

    // send B & recv merged B
    theChannel->sendVector(0, 0, *myVectB);
    theChannel->recvVector(0, 0, *vectB);
  }

  //
  // if main process, recv B & A from all, solve and send back X, B & result
  //

  else {

    *vectB = *myVectB;

    Vector remoteB(workArea, size);    
    // receive X and A contribution from subprocess & add them in

    for (int j=0; j<numChannels; j++) {

      Channel *theChannel = theChannels[j];
      theChannel->recvVector(0, 0, remoteB);
      *vectB += remoteB;
    }
  
    // send results back
    for (int j=0; j<numChannels; j++) {
      Channel *theChannel = theChannels[j];
      theChannel->sendVector(0, 0, *vectB);
    }
  }

#endif

  return *vectB;
}


int
MumpsParallelSOE::sendSelf(int commitTag, Channel &theChannel)
{
  int sendID =0;

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
	opserr << "MumpsParallelSOE::sendSelf() - failed to allocate channel array of size: " << 
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
	opserr << "MumpsParallelSOE::sendSelf() - failed to allocate id array of size: " << 
	  nextNumChannels << endln;
	return -1;
      }
      for (int i=0; i<numChannels; i++)
	localCol[i] = 0;    

      // allocate new processID for remote object
      sendID = numChannels;
    }

  } else 
    sendID = processID;


  // send remotes processID
  ID idData(2);
  idData(0) = sendID;
  idData(1) = matType;
  
  int res = theChannel.sendID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING MumpsParallelSOE::sendSelf() - failed to send data\n";
    return -1;
  }

  LinearSOESolver *theSoeSolver = this->getSolver();
  if (theSoeSolver != 0) {
    if (theSoeSolver->sendSelf(commitTag, theChannel) < 0) {
      opserr <<"WARNING MumpsParallelSOE::sendSelf() - failed to send solver\n";
      return -1;
    } 
  } else {
    opserr <<"WARNING MumpsParallelSOE::sendSelf() - no solver to send!\n";
    return -1;
  }
    

  return 0;
}


int 
MumpsParallelSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  ID idData(2);
  int res = theChannel.recvID(0, commitTag, idData);
  if (res < 0) {
    opserr <<"WARNING MumpsParallelSOE::recvSelf() - failed to send data\n";
    return -1;
  }	      
  processID = idData(0);
  matType = idData(1);

  numChannels = 1;
  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;

  localCol = new ID *[numChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;


  MumpsParallelSolver *theSolvr = new MumpsParallelSolver();
  if (theSolvr->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr <<"WARNING MumpsParallelSOE::sendSelf() - failed to recv solver\n";
    return -1;
  }
  
  theSolvr->setLinearSOE(*this);
  this->setSolver(*theSolvr);

  return 0;
}


int
MumpsParallelSOE::setProcessID(int dTag) 
{
  processID = dTag;
  return 0;
}

int
MumpsParallelSOE::setChannels(int nChannels, Channel **theC)
{
  numChannels = nChannels;

  if (theChannels != 0)
    delete [] theChannels;

  theChannels = new Channel *[numChannels];
  for (int i=0; i<numChannels; i++)
    theChannels[i] = theC[i];


  localCol = new ID *[nChannels];
  for (int i=0; i<numChannels; i++)
    localCol[i] = 0;

  return 0;
}
