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
                                                                        
// $Revision: 1.6 $
// $Date: 2009-05-11 20:56:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 02/06

// Description: This file contains the implementation for MumpsSOE

#include <MumpsSOE.h>
#include <MumpsSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>

#include <stdlib.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

// entries of an element matrix that addA() could not place in the sparsity
// pattern. Reported once when it first happens - a per-call message could be
// printed millions of times - and totalled when the SOE goes away.
static long long opsMumpsAddADropped = 0;
static bool      opsMumpsAddADroppedNoted = false;

MumpsSOE::MumpsSOE(MumpsSolver &the_Solver, int _matType)
:LinearSOE(the_Solver, LinSOE_TAGS_MumpsSOE),
 size(0), nnz(0), 
 A(0), B(0), X(0), 
 colA(0), rowA(0), rowB(0), colStartA(0),
 vectX(0), vectB(0),
 Asize(0), Bsize(0),
 factored(false), matType(_matType)
{
  the_Solver.setLinearSOE(*this);
}

MumpsSOE::MumpsSOE()
 :LinearSOE(LinSOE_TAGS_MumpsSOE),
  size(0), nnz(0), 
  A(0), B(0), X(0), 
  colA(0), rowA(0), rowB(0), colStartA(0),
  vectX(0), vectB(0),
  Asize(0), Bsize(0),
  factored(false), matType(0)
{

}

MumpsSOE::MumpsSOE(int classTag)
 :LinearSOE(classTag),
  size(0), nnz(0), 
  A(0), B(0), X(0), 
  colA(0), rowA(0), rowB(0), colStartA(0),
  vectX(0), vectB(0),
  Asize(0), Bsize(0),
  factored(false), matType(0)
{

}


MumpsSOE::MumpsSOE(LinearSOESolver &the_Solver, int classTag, int _matType)
  :LinearSOE(the_Solver, classTag),
   size(0), nnz(0), 
   A(0), B(0), X(0), 
   colA(0), rowA(0), rowB(0), colStartA(0),
   vectX(0), vectB(0),
   Asize(0), Bsize(0),
   factored(false), matType(_matType)
{

}


MumpsSOE::~MumpsSOE()
{
    if (opsMumpsAddADropped > 0) {
      opserr << "WARNING:MumpsSOE : addA() dropped " << (double)opsMumpsAddADropped
	     << " entr(ies) outside the sparsity pattern over the life of this"
	     << " system - the results it produced are not trustworthy\n";
      opsMumpsAddADropped = 0;
      opsMumpsAddADroppedNoted = false;
    }


    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (colStartA != 0) delete [] colStartA;
    if (rowA != 0) delete []rowA;
    if (colA != 0) delete []colA;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;
}


int
MumpsSOE::getNumEqn(void) const
{
    return size;
}

int 
MumpsSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int oldSize = size;
  size = theGraph.getNumVertex();
  
  // fist itearte through the vertices of the graph to get nnz
  Vertex *theVertex;
  int newNNZ = 0;
  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    const ID &theAdjacency = theVertex->getAdjacency();
    newNNZ += theAdjacency.Size() +1; // the +1 is for the diag entry
  }

  if (matType !=  0) {
    newNNZ -= size;
    newNNZ /= 2;
    newNNZ += size;
  }

  nnz = newNNZ;
  
  if (newNNZ > Asize) { // we have to get more space for A and rowA
    if (A != 0) delete [] A;
    if (rowA != 0) delete [] rowA;
    
    A = new double[newNNZ];
    rowA = new int[newNNZ];
    colA = new int[newNNZ];
    
    if (A == 0 || rowA == 0 || colA == 0) {
      opserr << "WARNING MumpsSOE::MumpsSOE :";
      opserr << " ran out of memory for A and rowA with nnz = ";
      opserr << newNNZ << " \n";
      size = 0; Asize = 0; nnz = 0;
      result =  -1;
    } 
    
    Asize = newNNZ;
  }
  
  // zero the matrix
  for (int i=0; i<Asize; i++)
    A[i] = 0;
  
  factored = false;
  
  if (size > Bsize) { // we have to get space for the vectors
    
    // delete the old	
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (colStartA != 0) delete [] colStartA;
    
    // create the new
    B = new double[size];
    X = new double[size];
    colStartA = new int[size+1]; 
    
    if (B == 0 || X == 0 || colStartA == 0) {
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
  }
  
  // create new Vectors objects
  if (size != oldSize) {
    if (vectX != 0)
      delete vectX;
    
    if (vectB != 0)
      delete vectB;
    
    vectX = new Vector(X,size);
    vectB = new Vector(B,size);	
  }
  
  // fill in colStartA and rowA
  if (size != 0) {
    colStartA[0] = 0;
    int startLoc = 0;
    int lastLoc = 0;
    for (int a=0; a<size; a++) {
      
      theVertex = theGraph.getVertexPtr(a);
      if (theVertex == 0) {
	opserr << "WARNING:MumpsSOE::setSize :";
	opserr << " vertex " << a << " not in graph! - size set to 0\n";
	size = 0;
	return -1;
      }
      
      int vertexTag = theVertex->getTag();
      rowA[lastLoc++] = vertexTag; // place diag in first
      const ID &theAdjacency = theVertex->getAdjacency();
      int idSize = theAdjacency.Size();
      
      // now we have to place the entries in the ID into order in rowA

      if (matType != 0) {

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

      colStartA[a+1] = lastLoc;;	    
      startLoc = lastLoc;
    }
  }

  // fill in colA
  int count = 0;
  for (int i=0; i<size; i++)
    for (int k=colStartA[i]; k<colStartA[i+1]; k++)
      colA[count++] = i;

  if (this->verifyRowAOrder() < 0)
    return -1;


  // invoke setSize() on the Solver    
  LinearSOESolver *the_Solver = this->getSolver();
  int solverOK = the_Solver->setSize();
  if (solverOK < 0) {
    opserr << "WARNING:MumpsSOE::setSize :";
    opserr << " solver failed setSize()\n";
    return solverOK;
  }    
  
  return result;
}

int
MumpsSOE::verifyRowAOrder(void)
{
  // addA() merges into the rows of a column, so the ascending order the setSize()
  // loops produce is a precondition, not a detail: the insertion sort places the
  // diagonal at the start of the column and every adjacency in its sorted slot,
  // in both matType branches and in the subclasses that build rowA themselves.
  // Checked rather than trusted - O(nnz) once per setSize against the O(L^2)
  // insertion sort it follows, and a violation would otherwise surface as a
  // silently wrong tangent instead of an error.
  for (int i=0; i<size; i++) {
    for (int k=colStartA[i]+1; k<colStartA[i+1]; k++) {
      if (rowA[k-1] >= rowA[k]) {
	opserr << "WARNING:MumpsSOE::verifyRowAOrder : rowA is not ascending in"
	       << " column " << i << " (rowA[" << k-1 << "]=" << rowA[k-1]
	       << " >= rowA[" << k << "]=" << rowA[k]
	       << ") - addA() cannot locate its entries\n";
	size = 0;
	return -1;
      }
    }
  }

  return 0;
}

int
MumpsSOE::addA(const Matrix &m, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  
	return 0;

    int idSize = id.Size();
    
    // check that m and id are of similar size
    if (idSize != m.noRows() && idSize != m.noCols()) {
	opserr << "MumpsSOE::addA() ";
	opserr << " - Matrix and ID not of similar sizes\n";
	return -1;
    }

    // The rows of a column are stored in ascending order (setSize() checks it),
    // so if the rows this element contributes are visited in ascending order too,
    // the slots of a whole column are located by ONE merge pass over the column
    // instead of one scan of it per row. The cost of a column drops from
    // O(idSize * L) to O(idSize + L), L being the column length: at L = 72, the
    // mean for a brick mesh, that is the difference between ~36 and ~1 probes per
    // entry. Bisection was the other candidate and was rejected: at these L it is
    // a dependent-load chain with unpredictable branches against a scan the
    // hardware prefetches perfectly, so it wins little or nothing.
    //
    // The visiting order is a permutation of the rows, not of the accumulations:
    // within one column each row appears once, so each A[k] still receives
    // exactly one addition per call, and the sort below is stable, so even a
    // repeated dof in `id` keeps the original relative order. The result is
    // therefore bit-identical to the scan it replaces.
    //
    // The dof numbers and the positions they came from are kept side by side so
    // the merge reads both with one load each instead of chasing id() through the
    // permutation.
    int ordRowStack[64], ordPosStack[64];
    int *ordRow = ordRowStack;
    int *ordPos = ordPosStack;
    int *ordHeap = 0;
    if (idSize > 64) {
      ordHeap = new int[2*idSize];
      ordRow = ordHeap;
      ordPos = ordHeap + idSize;
    }

    // stable insertion sort of the dofs of the element: idSize is the number of
    // dofs of one element, so this is a handful of entries, and its cost is paid
    // once per element and amortized over the idSize columns
    int nOrd = 0;
    for (int j=0; j<idSize; j++) {
      int row = id(j);
      if (row < 0 || row >= size)   // unnumbered or constrained: nothing to add
	continue;
      int p = nOrd++;
      while (p > 0 && ordRow[p-1] > row) {
	ordRow[p] = ordRow[p-1];
	ordPos[p] = ordPos[p-1];
	p--;
      }
      ordRow[p] = row;
      ordPos[p] = j;
    }

    const bool unit = (fact == 1.0);   // keep the arithmetic of fact == 1 exact

    for (int i=0; i<idSize; i++) {
      int col = id(i);
      if (col >= size || col < 0)
	continue;

      int k = colStartA[col];
      const int endColLoc = colStartA[col+1];

      for (int t=0; t<nOrd; t++) {
	const int row = ordRow[t];
	const int j   = ordPos[t];

	// matType != 0: only the triangle with row >= col is stored
	if (matType != 0 && row < col)
	  continue;

	// advance to the first stored row that is not before the wanted one. k is
	// never rewound: the targets arrive in ascending order, which is what
	// makes the whole column cost one pass. Equal rows leave k where it is, so
	// a repeated dof still finds its slot.
	while (k < endColLoc && rowA[k] < row)
	  k++;

	if (k < endColLoc && rowA[k] == row) {
	  if (unit)
	    A[k] += m(j,i);
	  else
	    A[k] += fact * m(j,i);
	} else {
	  // the entry is not in the sparsity pattern: it used to be dropped
	  // without a word, which turns into a wrong tangent and a wrong answer
	  // with nothing in the log to point at it
	  opsMumpsAddADropped++;
	  if (opsMumpsAddADroppedNoted == false) {
	    opsMumpsAddADroppedNoted = true;
	    opserr << "WARNING:MumpsSOE::addA : entry (row " << row << ", col "
		   << col << ") of an element matrix is not in the sparsity"
		   << " pattern - it is dropped, so the tangent is incomplete."
		   << " Further occurrences are counted, not printed.\n";
	  }
	}
      }  // for t
    }  // for i

    if (ordHeap != 0)
      delete [] ordHeap;

    return 0;
}


int
MumpsSOE::addB(const Vector &v, const ID &id, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;

    int idSize = id.Size();    
    // check that m and id are of similar size
    if (idSize != v.Size() ) {
	opserr << "MumpsSOE::addB() ";
	opserr << " - Vector and ID not of similar sizes\n";
	return -1;
    }    

    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<idSize; i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	

    return 0;
}


int
MumpsSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (v.Size() != size) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<size; i++) {
	    B[i] = v(i);
	}
    } else if (fact == -1.0) {
	for (int i=0; i<size; i++) {
	    B[i] = -v(i);
	}
    } else {
	for (int i=0; i<size; i++) {
	    B[i] = v(i) * fact;
	}
    }	

    return 0;
}

void 
MumpsSOE::zeroA(void)
{
    double *Aptr = A;
    for (int i=0; i<nnz; i++)
	*Aptr++ = 0;

	factored = false;
}
	
void 
MumpsSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}

void 
MumpsSOE::setX(int loc, double value)
{
    if (loc < size && loc >=0)
	X[loc] = value;
}

void 
MumpsSOE::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

const Vector &
MumpsSOE::getX(void)
{
    if (vectX == 0) {
	opserr << "FATAL MumpsSOE::getX - vectX == 0";
	exit(-1);
    }
    return *vectX;
}

const Vector &
MumpsSOE::getB(void)
{
    if (vectB == 0) {
	opserr << "FATAL MumpsSOE::getB - vectB == 0";
	exit(-1);
    }        
    return *vectB;
}

double 
MumpsSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
    
}    

int 
MumpsSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int 
MumpsSOE::recvSelf(int cTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)  
{
    return 0;
}

