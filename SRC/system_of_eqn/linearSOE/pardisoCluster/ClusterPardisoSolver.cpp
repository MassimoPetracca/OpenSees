/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: see ClusterPardisoSolver.h

#include <mpi.h>
#include <mkl_cluster_sparse_solver.h>

#include <ClusterPardisoSolver.h>
#include <ClusterPardisoSOE.h>
#include <classTags.h>
#include <OPS_Globals.h>
#include <stdlib.h>

ClusterPardisoSolver::ClusterPardisoSolver()
  : LinearSOESolver(SOLVER_TAGS_ClusterPardisoSolver),
    theClusterSOE(0), mtype(11), maxfct(1), mnum(1), msglvl(0),
    analyzed(false)
{
  for (int i = 0; i < 64; i++) {
    pt[i] = 0;
    iparm[i] = 0;
  }
  // defaults per the shipped MKL cluster examples (cl_solver_unsym_distr_c)
  iparm[0]  = 1;   // user-supplied iparm
  iparm[1]  = 2;   // METIS fill-in reordering
  iparm[7]  = 0;   // no iterative refinement (comparable to Mumps usage)
  iparm[9]  = 13;  // pivot perturbation 1e-13
  iparm[10] = 1;   // nonsymmetric scaling
  iparm[12] = 1;   // maximum weighted matching
  iparm[26] = (getenv("OPS_CPARDISO_DEBUG") != 0) ? 1 : 0;  // input check only in debug
  iparm[34] = 0;   // 1-based ia/ja
  iparm[39] = 2;   // distributed matrix, RHS and solution (row ranges)

  fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
}

ClusterPardisoSolver::~ClusterPardisoSolver()
{
  if (theClusterSOE != 0 && analyzed) {
    // release all internal memory
    int phase = -1;
    this->callSolver(phase, 0, 0);
  }
}

int
ClusterPardisoSolver::setLinearSOE(ClusterPardisoSOE &theSOE)
{
  theClusterSOE = &theSOE;
  return 0;
}

int
ClusterPardisoSolver::setIparm(int index, int value)
{
  if (index < 0 || index > 63) {
    opserr << "ClusterPardisoSolver::setIparm - index " << index
	   << " out of range [0,63] (0-based, MKL C documentation)\n";
    return -1;
  }
  // the integration owns these: 0 (user-iparm switch), 34 (0/1-based
  // indexing), 39-41 (distributed row-range layout). Overriding them
  // breaks the SOE<->solver contract in non-obvious ways.
  if (index == 0 || index == 34 || (index >= 39 && index <= 41)) {
    opserr << "ClusterPardisoSolver::setIparm - iparm[" << index
	   << "] is reserved by the OpenSees integration and cannot be set\n";
    return -1;
  }
  if (analyzed) {
    opserr << "ClusterPardisoSolver::setIparm - solver already analyzed, "
	   << "set options before the first analysis\n";
    return -1;
  }
  iparm[index] = value;
  return 0;
}

void
ClusterPardisoSolver::setMsglvl(int level)
{
  msglvl = level;
}

int
ClusterPardisoSolver::setMatrixType(int mt)
{
  if (mt != 11 && mt != 2 && mt != -2) {
    opserr << "ClusterPardisoSolver::setMatrixType - supported types: "
	   << "11 (unsymmetric), 2 (SPD), -2 (symmetric indefinite)\n";
    return -1;
  }
  if (analyzed) {
    opserr << "ClusterPardisoSolver::setMatrixType - solver already "
	   << "analyzed, set options before the first analysis\n";
    return -1;
  }
  mtype = mt;
  if (mt == 2 || mt == -2) {
    // scaling and weighted matching are unsymmetric-only features
    // (their defaults belong to mtype 11); leaving them on with a
    // symmetric type is at best ignored, at worst an input error
    iparm[10] = 0;
    iparm[12] = 0;
  }
  return 0;
}

// human-readable translation of the cluster_sparse_solver error codes
// (MKL Developer Reference); the raw code is printed alongside
static const char *
pardisoErrorText(int error)
{
  switch (error) {
  case -1:  return "inconsistent input (check matrix/iparm)";
  case -2:  return "not enough memory";
  case -3:  return "reordering problem";
  case -4:  return "zero pivot: numerical factorization or refinement problem (singular or near-singular matrix?)";
  case -5:  return "unclassified internal error";
  case -6:  return "reordering failed";
  case -7:  return "diagonal matrix is singular";
  case -8:  return "32-bit integer overflow (matrix too large for lp64 interface)";
  case -9:  return "not enough memory for out-of-core";
  case -10: return "problems opening out-of-core temporary files";
  case -11: return "read/write error on out-of-core files";
  default:  return "see MKL Developer Reference, cluster_sparse_solver error codes";
  }
}

int
ClusterPardisoSolver::callSolver(int phase, double *b, double *x)
{
  int n = theClusterSOE->size;
  int nrhs = 1;
  int idum = 0;
  int error = 0;
  double ddum = 0.0;

  // 1-based inclusive row range of this rank; empty domain: begin > end
  iparm[40] = theClusterSOE->rowBegin + 1;
  iparm[41] = theClusterSOE->rowEnd;

  double *aval = theClusterSOE->csrVal.empty() ? &ddum : theClusterSOE->csrVal.data();
  int *ia = theClusterSOE->csrIa.empty() ? &idum : theClusterSOE->csrIa.data();
  int *ja = theClusterSOE->csrJa.empty() ? &idum : theClusterSOE->csrJa.data();

  cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n,
			aval, ia, ja, &idum, &nrhs, iparm, &msglvl,
			(b != 0 ? (void *)b : (void *)&ddum),
			(x != 0 ? (void *)x : (void *)&ddum),
			&fcomm, &error);

  if (error != 0)
    opserr << "ERROR ClusterPardisoSolver phase " << phase
	   << " returned error " << error << ": "
	   << pardisoErrorText(error) << "\n";
  return error;
}

int
ClusterPardisoSolver::setSize(void)
{
  if (theClusterSOE == 0) {
    opserr << "ClusterPardisoSolver::setSize - no associated SOE\n";
    return -1;
  }

  if (analyzed) {
    // new structure: release the previous internal storage first
    int phase = -1;
    this->callSolver(phase, 0, 0);
    analyzed = false;
  }

  // phase 11 is DEFERRED to the first solve(): at setSize the matrix
  // values are still zero, and with matching/scaling on (iparm[10],
  // iparm[12]) the analysis uses the numerical values -- run on a zero
  // matrix it produces garbage scale factors and a factorization of the
  // wrong (mis-scaled) matrix. MUMPS does not have this trap: its job=1
  // is purely structural.
  return 0;
}

int
ClusterPardisoSolver::solve(void)
{
  if (theClusterSOE == 0) {
    opserr << "ClusterPardisoSolver::solve - no associated SOE\n";
    return -1;
  }
  int error = 0;

  if (!analyzed) {
    // deferred analysis: first solve after a (re)setSize, values are real
    theClusterSOE->packOwnedValues();
    error = this->callSolver(11, 0, 0);
    if (error != 0)
      return -2;
    analyzed = true;
  }

  if (theClusterSOE->factored == false) {
    theClusterSOE->packOwnedValues();
    error = this->callSolver(22, 0, 0);
    if (error != 0)
      return -3;
    theClusterSOE->factored = true;
  }

  // distributed solve: b and x are the OWNED slices of the full mirrors
  double *bSlice = theClusterSOE->B + theClusterSOE->rowBegin;
  double *xSlice = theClusterSOE->X + theClusterSOE->rowBegin;
  error = this->callSolver(33, bSlice, xSlice);
  if (error != 0)
    return -4;
  return 0;
}

int
ClusterPardisoSolver::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int
ClusterPardisoSolver::recvSelf(int cTag, Channel &theChannel,
			       FEM_ObjectBroker &theBroker)
{
  return 0;
}
