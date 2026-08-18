/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: LinearSOESolver wrapping MKL cluster_sparse_solver
// (distributed assembled CSR input, iparm[39]=2). Phase mapping:
//   11 (analysis/reordering)  once per structure, at the FIRST solve()
//      after setSize -- NOT in setSize itself: with matching/scaling on
//      (iparm[10]/iparm[12]) the analysis uses the numerical values,
//      which are still zero at setSize time
//   22 (numerical factorization) at each new tangent (factored==false)
//   33 (solve)                per RHS
// mtype=11 (real unsymmetric) to match the OpenSees Mumps usage
// (matType 0); ia/ja and the row range are 1-based as in the MKL
// examples; communicator is MPI_COMM_WORLD via MPI_Comm_c2f.

#ifndef ClusterPardisoSolver_h
#define ClusterPardisoSolver_h

#include <LinearSOESolver.h>

class ClusterPardisoSOE;

class ClusterPardisoSolver : public LinearSOESolver
{
 public:
  ClusterPardisoSolver();
  ~ClusterPardisoSolver();

  int solve(void);
  int setSize(void);
  int setLinearSOE(ClusterPardisoSOE &theSOE);

  // user overrides from the Tcl flags (applied over the constructor
  // defaults, before the first cluster_sparse_solver call). setIparm
  // takes the 0-based C index of the MKL documentation and REFUSES the
  // indices the integration owns (0, 34, 39, 40, 41).
  int setIparm(int index, int value);
  void setMsglvl(int level);

  // matrix type: 11 (default, real unsymmetric LU), 2 (SPD, Cholesky),
  // -2 (symmetric indefinite, LDL^T). For +-2 the SOE hands over only
  // the upper triangle and scaling/matching are turned off (they are
  // unsymmetric-only features whose analysis reads the values).
  int setMatrixType(int mt);
  int getMatrixType(void) const { return mtype; }

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

 private:
  int callSolver(int phase, double *b, double *x);

  ClusterPardisoSOE *theClusterSOE;

  void *pt[64];
  int iparm[64];
  int mtype, maxfct, mnum, msglvl;
  int fcomm;       // fortran communicator handle
  bool analyzed;   // phase 11 done for the current structure
};

#endif
