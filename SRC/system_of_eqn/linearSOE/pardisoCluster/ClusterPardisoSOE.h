/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: ClusterPardisoSOE is a LinearSOE for OpenSeesMP built on
// MKL cluster_sparse_solver with distributed assembled input
// (iparm[39]=2): each rank owns the contiguous equation block reported
// by EquationPartition (ParallelNumberer Plain ordering) and provides
// its rows of A in row-major CSR (1-based, as in the MKL examples).
// Element contributions to rows owned by OTHER ranks are shipped to the
// owner (pattern once at setSize, values at each factorization).
// B and X keep full-length mirrors on every rank (MPI_Allreduce /
// MPI_Allgatherv) for compatibility with the replicated-vector
// ArpackSolver and the convergence tests; the storage/assembly scheme
// (COO grouped by column, from the local graph) is inherited from
// MumpsSOE unchanged.
//
// Cross-rank exchange: the pattern of the foreign entries is shipped to
// the owner once per setSize (MPI_Alltoallv of (row,col) pairs, send
// order sorted by (owner,row,col)); the merged owned CSR is the union
// of local-owned and received pairs, with scatter-add maps for both.
// Values travel at each factorization through a fixed-size Alltoallv in
// packOwnedValues -- which is therefore COLLECTIVE when np > 1 (always
// reached in lockstep: the solver phases are collective).

#ifndef ClusterPardisoSOE_h
#define ClusterPardisoSOE_h

#include <mpi.h>
#include <MumpsSOE.h>
#include <Vector.h>
#include <vector>

class ClusterPardisoSolver;

class ClusterPardisoSOE : public MumpsSOE
{
 public:
  ClusterPardisoSOE(ClusterPardisoSolver &theSolver);
  ~ClusterPardisoSOE();

  int setSize(Graph &theGraph);
  int solve(void);

  // local accumulation in myB; the assembled global B lives in the
  // base-class mirror B after solve() (same contract as MumpsParallelSOE)
  int addB(const Vector &, const ID &, double fact = 1.0);
  int setB(const Vector &, double fact = 1.0);
  const Vector &getB(void);
  void zeroB(void);

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  friend class ClusterPardisoSolver;

 private:
  int buildOwnedCSR(void);  // split local COO into owned/foreign, build CSR + maps
  void packOwnedValues(void);

  int rank, np;

  // typed pointer to the solver: buildOwnedCSR reads getMatrixType()
  // to decide full vs upper-triangle CSR (symmetric mtypes)
  ClusterPardisoSolver *thePardisoSolver;

  double *myB;
  Vector *myVectB;

  // owned global row range [rowBegin, rowEnd) and per-rank layout for the
  // X allgather (counts/displs in rank order)
  int rowBegin, rowEnd;
  std::vector<int> gatherCounts, gatherDispls;

  // owned-row CSR (1-based ia/ja, row-major) handed to cluster_sparse_solver:
  // union of the local owned entries and the pairs received from other
  // ranks, one slot per distinct (row,col)
  std::vector<int> csrIa, csrJa;
  std::vector<double> csrVal;
  std::vector<int> ownedCOO;  // A[] indices of the locally-owned COO entries
  std::vector<int> cooToCsr;  // CSR value slot of the m-th ownedCOO entry

  // cross-rank exchange (fixed after setSize):
  // foreignCOO: A[] indices of entries whose row belongs to another rank,
  // sorted by (owner,row,col) = the value-send order; recvToCsr: CSR slot
  // of each received value; counts/displs in doubles for the Alltoallv
  std::vector<int> foreignCOO;
  std::vector<int> recvToCsr;
  std::vector<int> valSendCounts, valSendDispls, valRecvCounts, valRecvDispls;
  std::vector<double> sendVals, recvVals;
};

#endif
