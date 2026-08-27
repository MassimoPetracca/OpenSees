/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Written: 2026/08
//
// The collective agreement primitives of OpenSeesMP. See ParallelAgreement.h for
// what each one promises and, more importantly, for why the predicate is not
// `np > 1`.

#include <ParallelAgreement.h>
#include <OPS_Globals.h>

#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif // _PARALLEL_INTERPRETERS

// Raised by the objects only a coupled run can have - see the header. Not lowered.
static bool ops_coupled_parallel_run = false;

// How many agreement points this process has passed. Folded into the reduce of
// OPS_agreeWorst, so it costs no collective of its own.
static int ops_agreement_point = 0;

// WHO ORIGINATED THE FAILURE, which stops being knowable from the return codes the
// moment OPS_agreeWorst() starts handing every process the same one. Both describe
// the phase currently running and are cleared by OPS_agreeWorstWho(), which ends it.
static bool ops_failed_here = false;      // this process produced a negative code
static bool ops_failed_elsewhere = false; // it only learned of one

void
OPS_setCoupledParallelRun(void)
{
    ops_coupled_parallel_run = true;
}

bool
OPS_inCoupledParallelRun(void)
{
#ifdef _PARALLEL_INTERPRETERS
    if (ops_coupled_parallel_run) {
      int np = 1;
      MPI_Comm_size(MPI_COMM_WORLD, &np);
      return (np > 1);
    }
#endif // _PARALLEL_INTERPRETERS

    return false;
}

int
OPS_agreementRank(void)
{
#ifdef _PARALLEL_INTERPRETERS
    if (OPS_inCoupledParallelRun()) {
      int here = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &here);
      return here;
    }
#endif // _PARALLEL_INTERPRETERS

    return 0;
}

#ifdef _PARALLEL_INTERPRETERS
// The processes are no longer executing the same sequence of agreement points, so
// nothing after this can be trusted: every subsequent collective pairs a call of
// one process with a different call of another. Reported by every process, because
// the useful part is each one's own count.
//
// This fails the analysis rather than calling MPI_Abort - nothing in SRC aborts,
// and a script that gets a failure back can at least close its recorders. The run
// is over either way.
static int
ops_lostLockstep(int seqHere, int seqMin, int seqMax)
{
    int here = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &here);

    opserr << "FATAL - the processes are no longer in step: process " << here
	   << " is at agreement point " << seqHere
	   << " while the others are between " << seqMin << " and " << seqMax << ".\n"
	   << "        An earlier phase returned on one process without agreeing with "
	   << "the others.\n        The analysis cannot continue.\n";

    return OPS_AGREEMENT_LOST;
}
#endif // _PARALLEL_INTERPRETERS

int
OPS_agreeWorst(int resultHere)
{
#ifdef _PARALLEL_INTERPRETERS
    if (OPS_inCoupledParallelRun()) {

      // +seq and -seq under one MPI_MIN give both the smallest and the largest
      // count over the processes, so the lockstep check rides along free
      const int seq = ops_agreement_point++;

      if (resultHere < 0)
	ops_failed_here = true;

      int in[3], out[3];
      in[0] = resultHere;  in[1] = seq;  in[2] = -seq;
      out[0] = resultHere; out[1] = seq; out[2] = -seq;

      if (MPI_Allreduce(in, out, 3, MPI_INT, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS) {
	opserr << "OPS_agreeWorst - MPI_Allreduce failed\n";
	return resultHere;
      }

      if (out[1] != -out[2])
	return ops_lostLockstep(seq, out[1], -out[2]);

      if (out[0] < 0 && resultHere >= 0)
	ops_failed_elsewhere = true;

      return out[0];
    }
#endif // _PARALLEL_INTERPRETERS

    return resultHere;
}

int
OPS_agreeWorstWho(int resultHere, int &worstRank, int &numFailed, int &numProcesses)
{
#ifdef _PARALLEL_INTERPRETERS
    if (OPS_inCoupledParallelRun()) {

      int np = 1, here = 0;
      MPI_Comm_size(MPI_COMM_WORLD, &np);
      MPI_Comm_rank(MPI_COMM_WORLD, &here);

      ops_agreement_point++;

      int in[2], out[2];
      in[0] = resultHere;  in[1] = here;
      out[0] = resultHere; out[1] = here;

      if (MPI_Allreduce(in, out, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD) != MPI_SUCCESS) {
	opserr << "OPS_agreeWorstWho - MPI_Allreduce failed\n";
	return resultHere;
      }

      worstRank = out[1];
      numProcesses = np;
      numFailed = 0;

      if (out[0] < 0) {
	// Everything below is guarded on out[0], which is identical on every
	// process, so it is taken by all of them or by none. Guarding on
	// resultHere would not be.
	//
	// THE VERDICT IS ABOVE, THE ATTRIBUTION IS HERE, and they cannot be the
	// same reduce. By the time a phase ends, OPS_agreeWorst() has already given
	// every process the same negative code, so the MINLOC over the code alone
	// names the LOWEST-NUMBERED process rather than the one whose element
	// failed - it would have blamed process 0 for a material on process 1.
	// So the attribution reduces the code only where this process ORIGINATED
	// the failure, and 1 - never a failure code, so MINLOC cannot pick it -
	// where it is merely carrying someone else's verdict.
	int originatedHere = (ops_failed_here ||
			      (resultHere < 0 && !ops_failed_elsewhere)) ? 1 : 0;

	int who[2], whoOut[2];
	who[0] = originatedHere ? resultHere : 1;  who[1] = here;
	whoOut[0] = who[0];                        whoOut[1] = who[1];
	if (MPI_Allreduce(who, whoOut, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD)
	    == MPI_SUCCESS && whoOut[0] < 0)
	  worstRank = whoOut[1];

	// and how many originated it, which is what separates one partition's
	// element or material from a global problem, typically the solver
	if (MPI_Allreduce(&originatedHere, &numFailed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD)
	    != MPI_SUCCESS)
	  numFailed = 0;
      }

      // the phase is over, so what happened inside it stops being current
      ops_failed_here = false;
      ops_failed_elsewhere = false;

      return out[0];
    }
#endif // _PARALLEL_INTERPRETERS

    return resultHere;
}

bool
OPS_agreeAny(bool flagHere)
{
#ifdef _PARALLEL_INTERPRETERS
    if (OPS_inCoupledParallelRun()) {

      ops_agreement_point++;

      int here = flagHere ? 1 : 0;
      int anywhere = here;
      if (MPI_Allreduce(&here, &anywhere, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS) {
	opserr << "OPS_agreeAny - MPI_Allreduce failed\n";
	return flagHere;
      }

      return (anywhere != 0);
    }
#endif // _PARALLEL_INTERPRETERS

    return flagHere;
}
