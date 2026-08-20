/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: see ClusterPardisoSOE.h

#include <ClusterPardisoSOE.h>
#include <ClusterPardisoSolver.h>
#include <EquationPartition.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <ID.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <algorithm>
#include <array>
#include <stdlib.h>
#include <math.h>

ClusterPardisoSOE::ClusterPardisoSOE(ClusterPardisoSolver &theSolvr)
  : MumpsSOE(theSolvr, LinSOE_TAGS_ClusterPardisoSOE, 0),
    rank(0), np(1), thePardisoSolver(&theSolvr),
    myB(0), myVectB(0), rowBegin(0), rowEnd(0)
{
  theSolvr.setLinearSOE(*this);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
}

ClusterPardisoSOE::~ClusterPardisoSOE()
{
  if (myB != 0)
    delete [] myB;
  if (myVectB != 0)
    delete myVectB;
}

int
ClusterPardisoSOE::setSize(Graph &theGraph)
{
  int result = 0;
  int oldSize = size;

  // local sparsity: count nnz of the local graph (diag + adjacency),
  // exactly as MumpsParallelSOE does (matType is fixed to 0: full
  // unsymmetric pattern, PARDISO mtype 11)
  int maxVertexTag = -1;
  Vertex *theVertex;
  int newNNZ = 0;
  size = theGraph.getNumVertex();

  VertexIter &theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    int vertexTag = theVertex->getTag();
    if (vertexTag > maxVertexTag)
      maxVertexTag = vertexTag;
    const ID &theAdjacency = theVertex->getAdjacency();
    newNNZ += theAdjacency.Size() + 1;
  }
  nnz = newNNZ;

  // global size = allreduce-max of the local max equation number, +1
  int globalMaxTag = maxVertexTag;
  MPI_Allreduce(&maxVertexTag, &globalMaxTag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  size = globalMaxTag + 1;

  if (nnz > Asize) {
    if (A != 0) delete [] A;
    if (rowA != 0) delete [] rowA;
    if (colA != 0) delete [] colA;

    A = new double[nnz];
    rowA = new int[nnz];
    colA = new int[nnz];
    for (int i = 0; i < nnz; i++) {
      A[i] = 0.0; rowA[i] = 0; colA[i] = 0;
    }
    Asize = nnz;
  }

  if (size > Bsize) {
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (myB != 0) delete [] myB;
    if (colStartA != 0) delete [] colStartA;

    B = new double[size];
    X = new double[size];
    myB = new double[size];
    colStartA = new int[size + 1];
    Bsize = size;
  }
  for (int j = 0; j < size; j++) {
    B[j] = 0.0; X[j] = 0.0; myB[j] = 0.0;
  }

  if (size != oldSize) {
    if (vectX != 0) delete vectX;
    if (vectB != 0) delete vectB;
    if (myVectB != 0) delete myVectB;
    vectX = new Vector(X, size);
    vectB = new Vector(B, size);
    myVectB = new Vector(myB, size);
  }

  // fill colStartA / rowA / colA: for each local column (vertex) the
  // sorted row indices (diag + adjacency) -- same layout as MumpsSOE,
  // global 0-based indices, columns without a local vertex stay empty
  if (size != 0) {
    colStartA[0] = 0;
    int startLoc = 0;
    int lastLoc = 0;
    for (int a = 0; a < size; a++) {
      theVertex = theGraph.getVertexPtr(a);
      if (theVertex != 0) {
	int vertexTag = theVertex->getTag();
	rowA[lastLoc++] = vertexTag;
	const ID &theAdjacency = theVertex->getAdjacency();
	int idSize = theAdjacency.Size();
	for (int i = 0; i < idSize; i++) {
	  int row = theAdjacency(i);
	  bool foundPlace = false;
	  for (int j = startLoc; j < lastLoc; j++)
	    if (rowA[j] > row) {
	      for (int k = lastLoc; k > j; k--)
		rowA[k] = rowA[k - 1];
	      rowA[j] = row;
	      foundPlace = true;
	      j = lastLoc;
	    }
	  if (foundPlace == false)
	    rowA[lastLoc] = row;
	  lastLoc++;
	}
      }
      colStartA[a + 1] = lastLoc;
      startLoc = lastLoc;
    }
    int count = 0;
    for (int i = 0; i < size; i++)
      for (int k = colStartA[i]; k < colStartA[i + 1]; k++)
	colA[count++] = i;
  }

  // ownership: EquationPartition from the ParallelNumberer Plain ordering;
  // at np==1 the whole range is trivially owned
  EquationPartition &part = EquationPartition::instance();
  if (np == 1) {
    rowBegin = 0;
    rowEnd = size;
  } else {
    if (!part.isValid() || part.getNumRanks() != np ||
	part.getNumEquations() != size) {
      opserr << "ERROR ClusterPardisoSOE::setSize - no valid equation "
	     << "partition; system ClusterPardiso requires numberer "
	     << "ParallelPlain in OpenSeesMP\n";
      return -1;
    }
    rowBegin = part.getMyRowBegin();
    rowEnd = part.getMyRowEnd();
  }

  // allgather layout for X (counts/displs by rank)
  gatherCounts.assign(np, 0);
  gatherDispls.assign(np, 0);
  for (int r = 0; r < np; r++) {
    int b = (np == 1) ? 0 : part.getRowBegin(r);
    int e = (np == 1) ? size : part.getRowEnd(r);
    gatherCounts[r] = e - b;
    gatherDispls[r] = b;
  }

  // buildOwnedCSR contains collectives (pattern exchange + status check):
  // every rank reaches it unconditionally, and its verdict is collective
  result = this->buildOwnedCSR();
  if (result < 0)
    return result;

  factored = false;

  LinearSOESolver *theSolvr = this->getSolver();
  int solverOK = theSolvr->setSize();
  if (solverOK < 0) {
    opserr << "WARNING ClusterPardisoSOE::setSize - solver failed setSize()\n";
    return solverOK;
  }
  return result;
}

int
ClusterPardisoSOE::buildOwnedCSR(void)
{
  EquationPartition &part = EquationPartition::instance();

  // symmetric mtypes (2, -2): MKL wants ONLY the upper triangle in the
  // CSR, diagonal included. Assembly stays full-pattern (its cost is
  // negligible next to the factorization); the lower triangle is simply
  // not handed over -- and not shipped across ranks either.
  bool upperOnly = (thePardisoSolver->getMatrixType() != 11);

  ownedCOO.clear();
  foreignCOO.clear();

  for (int k = 0; k < nnz; k++) {
    int row = rowA[k];
    if (upperOnly && colA[k] < row)
      continue;
    if (row >= rowBegin && row < rowEnd)
      ownedCOO.push_back(k);
    else
      foreignCOO.push_back(k);
  }

  // sort owned triplets by (row, col) -> row-major CSR
  std::sort(ownedCOO.begin(), ownedCOO.end(),
	    [this](int a, int b) {
	      if (rowA[a] != rowA[b]) return rowA[a] < rowA[b];
	      return colA[a] < colA[b];
	    });

  // ---- pattern exchange: ship each foreign (row,col) to the row's owner.
  // Send order = foreignCOO sorted by (owner,row,col): grouped by
  // destination rank and deterministic; the same order carries the VALUES
  // at every factorization (packOwnedValues), so it is frozen here.
  std::vector<int> recvPairs;
  if (np > 1) {
    std::sort(foreignCOO.begin(), foreignCOO.end(),
	      [this, &part](int a, int b) {
		int oa = part.getOwner(rowA[a]);
		int ob = part.getOwner(rowA[b]);
		if (oa != ob) return oa < ob;
		if (rowA[a] != rowA[b]) return rowA[a] < rowA[b];
		return colA[a] < colA[b];
	      });

    valSendCounts.assign(np, 0);
    for (size_t m = 0; m < foreignCOO.size(); m++)
      valSendCounts[part.getOwner(rowA[foreignCOO[m]])]++;

    valRecvCounts.assign(np, 0);
    MPI_Alltoall(valSendCounts.data(), 1, MPI_INT,
		 valRecvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    valSendDispls.assign(np, 0);
    valRecvDispls.assign(np, 0);
    std::vector<int> pairSendCounts(np), pairSendDispls(np);
    std::vector<int> pairRecvCounts(np), pairRecvDispls(np);
    int sTot = 0, rTot = 0;
    for (int r = 0; r < np; r++) {
      valSendDispls[r] = sTot;
      valRecvDispls[r] = rTot;
      pairSendCounts[r] = 2 * valSendCounts[r];
      pairRecvCounts[r] = 2 * valRecvCounts[r];
      pairSendDispls[r] = 2 * sTot;
      pairRecvDispls[r] = 2 * rTot;
      sTot += valSendCounts[r];
      rTot += valRecvCounts[r];
    }

    std::vector<int> sendPairs(2 * (size_t)sTot);
    for (size_t m = 0; m < foreignCOO.size(); m++) {
      sendPairs[2 * m]     = rowA[foreignCOO[m]];
      sendPairs[2 * m + 1] = colA[foreignCOO[m]];
    }
    recvPairs.resize(2 * (size_t)rTot);
    MPI_Alltoallv(sendPairs.data(), pairSendCounts.data(),
		  pairSendDispls.data(), MPI_INT,
		  recvPairs.data(), pairRecvCounts.data(),
		  pairRecvDispls.data(), MPI_INT, MPI_COMM_WORLD);

    sendVals.assign(sTot, 0.0);
    recvVals.assign(rTot, 0.0);
  } else {
    valSendCounts.clear(); valSendDispls.clear();
    valRecvCounts.clear(); valRecvDispls.clear();
    sendVals.clear(); recvVals.clear();
  }
  int nRecv = (int)recvPairs.size() / 2;

  // ---- merged owned pattern: union of local owned entries and received
  // pairs, one CSR slot per distinct (row,col). Each entry carries its
  // source so the scatter-add maps come out of the same walk:
  // src >= 0 -> ownedCOO index, src < 0 -> received slot -(src+1)
  std::vector<std::array<int, 3>> ents;
  ents.reserve(ownedCOO.size() + (size_t)nRecv);
  for (size_t m = 0; m < ownedCOO.size(); m++) {
    std::array<int, 3> t = { rowA[ownedCOO[m]], colA[ownedCOO[m]], (int)m };
    ents.push_back(t);
  }
  for (int i = 0; i < nRecv; i++) {
    std::array<int, 3> t = { recvPairs[2 * i], recvPairs[2 * i + 1], -i - 1 };
    ents.push_back(t);
  }
  std::sort(ents.begin(), ents.end());  // lexicographic: (row, col, src)

  int nOwnedRows = rowEnd - rowBegin;
  csrIa.assign(nOwnedRows + 1, 0);
  csrJa.clear();
  csrJa.reserve(ents.size());
  cooToCsr.assign(ownedCOO.size(), -1);
  recvToCsr.assign(nRecv, -1);

  int localFail = 0;
  csrIa[0] = 1;  // 1-based, MKL example convention
  size_t e = 0;
  for (int r = 0; r < nOwnedRows; r++) {
    int globalRow = rowBegin + r;
    int lastCol = -1;
    while (e < ents.size() && ents[e][0] == globalRow) {
      int col = ents[e][1];
      if (col != lastCol) {
	csrJa.push_back(col + 1);
	lastCol = col;
      }
      int src = ents[e][2];
      if (src >= 0)
	cooToCsr[src] = (int)csrJa.size() - 1;
      else
	recvToCsr[-src - 1] = (int)csrJa.size() - 1;
      e++;
    }
    csrIa[r + 1] = (int)csrJa.size() + 1;
  }
  if (e != ents.size()) {
    // a received (or local) row fell outside [rowBegin, rowEnd): partition
    // inconsistency between sender and receiver
    opserr << "ERROR ClusterPardisoSOE::buildOwnedCSR - rank " << rank
	   << ": " << (int)(ents.size() - e)
	   << " entries outside the owned row range\n";
    localFail = 1;
  }
  csrVal.assign(csrJa.size(), 0.0);

  // collective verdict: a rank that returns early while the others reach
  // the collectives in packOwnedValues/solve deadlocks the job
  int globalFail = 0;
  MPI_Allreduce(&localFail, &globalFail, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (globalFail > 0)
    return -3;
  return 0;
}

void
ClusterPardisoSOE::packOwnedValues(void)
{
  std::fill(csrVal.begin(), csrVal.end(), 0.0);
  for (size_t m = 0; m < ownedCOO.size(); m++)
    csrVal[cooToCsr[m]] += A[ownedCOO[m]];

  if (np > 1) {
    // value exchange, fixed layout from buildOwnedCSR (COLLECTIVE)
    for (size_t m = 0; m < foreignCOO.size(); m++)
      sendVals[m] = A[foreignCOO[m]];
    MPI_Alltoallv(sendVals.data(), valSendCounts.data(),
		  valSendDispls.data(), MPI_DOUBLE,
		  recvVals.data(), valRecvCounts.data(),
		  valRecvDispls.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    for (size_t i = 0; i < recvVals.size(); i++)
      csrVal[recvToCsr[i]] += recvVals[i];
  }
}

int
ClusterPardisoSOE::solve(void)
{
  // assembled global B on every rank (mirror), summing the per-rank myB
  MPI_Allreduce(myB, B, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  int resSolver = this->LinearSOE::solve();

  if (resSolver == 0) {
    // solver wrote the owned slice of X: complete the mirror
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
		   X, gatherCounts.data(), gatherDispls.data(),
		   MPI_DOUBLE, MPI_COMM_WORLD);
  }

  // temporary layer-separation diagnostics: residual of the returned x
  // against (1) the assembled local COO and (2) the CSR handed to PARDISO
  if (getenv("OPS_CPARDISO_DEBUG") != 0) {
    double rnCOO = 0.0, rnCSR = 0.0, bn = 0.0, xn = 0.0;
    std::vector<double> r(size, 0.0);
    for (int k = 0; k < nnz; k++)
      r[rowA[k]] += A[k] * X[colA[k]];
    // the local COO holds only this rank's element contributions: the sum
    // across ranks is the global A*x (B is already the assembled mirror)
    MPI_Allreduce(MPI_IN_PLACE, r.data(), size, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);
    for (int i = 0; i < size; i++) {
      rnCOO += (r[i] - B[i]) * (r[i] - B[i]);
      bn += B[i] * B[i];
      xn += X[i] * X[i];
    }
    // the CSR check is only meaningful for the full pattern: with a
    // symmetric mtype the CSR holds the upper triangle alone and the
    // row sums are incomplete by construction
    bool upperOnly = (thePardisoSolver->getMatrixType() != 11);
    int nOwnedRows = rowEnd - rowBegin;
    for (int rr = 0; rr < nOwnedRows && !upperOnly; rr++) {
      double ax = 0.0;
      for (int s = csrIa[rr] - 1; s < csrIa[rr + 1] - 1; s++)
	ax += csrVal[s] * X[csrJa[s] - 1];
      double d = ax - B[rowBegin + rr];
      rnCSR += d * d;
    }
    opserr << "CPARDISO DEBUG rank " << rank << ": n " << size
	   << " nnz " << nnz << " csrNnz " << (int)csrJa.size()
	   << " ||Acoo*x-b|| " << sqrt(rnCOO);
    if (upperOnly)
      opserr << " ||Acsr*x-b|| n/a(sym)";
    else
      opserr << " ||Acsr*x-b||(owned) " << sqrt(rnCSR);
    opserr << " ||b|| " << sqrt(bn) << " ||x|| " << sqrt(xn) << "\n";
  }
  return resSolver;
}

int
ClusterPardisoSOE::addB(const Vector &v, const ID &id, double fact)
{
  if (fact == 0.0)
    return 0;
  int idSize = id.Size();
  if (idSize != v.Size()) {
    opserr << "ClusterPardisoSOE::addB - Vector and ID not of similar sizes\n";
    return -1;
  }
  if (fact == 1.0) {
    for (int i = 0; i < idSize; i++) {
      int pos = id(i);
      if (pos < size && pos >= 0)
	myB[pos] += v(i);
    }
  } else if (fact == -1.0) {
    for (int i = 0; i < idSize; i++) {
      int pos = id(i);
      if (pos < size && pos >= 0)
	myB[pos] -= v(i);
    }
  } else {
    for (int i = 0; i < idSize; i++) {
      int pos = id(i);
      if (pos < size && pos >= 0)
	myB[pos] += v(i) * fact;
    }
  }
  return 0;
}

int
ClusterPardisoSOE::setB(const Vector &v, double fact)
{
  if (fact == 0.0)
    return 0;
  if (v.Size() != size) {
    opserr << "WARNING ClusterPardisoSOE::setB - incompatible sizes "
	   << size << " and " << v.Size() << endln;
    return -1;
  }
  if (fact == 1.0) {
    for (int i = 0; i < size; i++)
      myB[i] = v(i);
  } else {
    for (int i = 0; i < size; i++)
      myB[i] = v(i) * fact;
  }
  return 0;
}

const Vector &
ClusterPardisoSOE::getB(void)
{
  // refresh the assembled global mirror from the per-rank myB, so callers
  // that read B WITHOUT a preceding solve() get the true assembled RHS.
  // DistributedDisplacementControl::domainChanged() probes getB() to build
  // its reference-load vector phat before any solve; without this refresh B
  // is still the setSize() zero -> "zero reference load" and the analysis
  // aborts. Collective on MPI_COMM_WORLD (size is the global equation count,
  // identical on every rank), same collective contract as
  // MumpsParallelSOE::getB(); mirrors the Allreduce at the top of solve().
  if (size > 0)
    MPI_Allreduce(myB, B, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return *vectB;
}

void
ClusterPardisoSOE::zeroB(void)
{
  for (int i = 0; i < size; i++)
    myB[i] = 0.0;
}

int
ClusterPardisoSOE::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int
ClusterPardisoSOE::recvSelf(int cTag, Channel &theChannel,
			    FEM_ObjectBroker &theBroker)
{
  return 0;
}
