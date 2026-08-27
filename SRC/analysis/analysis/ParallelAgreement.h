/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

#ifndef ParallelAgreement_h
#define ParallelAgreement_h

// Written: 2026/08
//
// THE COLLECTIVE AGREEMENT PRIMITIVES OF A PARALLEL INTERPRETER RUN (OpenSeesMP),
// and the one predicate that says whether a collective is owed at all.
//
// They live in a file of their own so that every process-agreement decision in the
// analysis is made by the same three functions and gated by the same predicate.
// The alternative - each caller writing its own np/rank/flag logic - is what left
// rank 0 out of an Allreduce that every other rank had entered
// (see CTestImplexWrapper::aggregateImplexError).
//
// WHY A PREDICATE AND NOT `np > 1`. OpenSeesMP is used in two completely different
// ways, and only one of them owes anyone a collective:
//
//   * ONE model partitioned over the processes. The equations are numbered
//     globally, solve() and getB() are collectives, and a process that abandons a
//     step alone hangs the job. This is the case the agreement exists for.
//   * N INDEPENDENT ANALYSES, one per process - parametric studies, the other
//     reason people reach for OpenSeesMP. The processes share nothing, run
//     different numbers of steps, and an unconditional collective there is worse
//     than a hang: MPI pairs the reduce of one process's step with the reduce of
//     ANOTHER process's unrelated step, so one analysis's failure silently
//     abandons a step of a different analysis.
//
// The discriminator is exact and costs nothing. A coupled run MUST have a
// ParallelNumberer - without one every process numbers its own equations from 0
// and there is no global system to solve - so the flag is raised by the objects
// that only a coupled run can have: ParallelNumberer and the parallel SOEs, in
// their setChannels(). No front end has to be told anything.
//
// It is never lowered. A script that ran a coupled analysis and then switched to
// independent ones in the SAME process would keep paying the collectives; if that
// ever matters, the place to lower it is the `wipe` command.
bool OPS_inCoupledParallelRun(void);

// Raises the flag above. Called by ParallelNumberer::setChannels() and by the
// parallel SOEs. Idempotent.
void OPS_setCoupledParallelRun(void);

// What OPS_agreeWorst returns when the processes have lost lockstep: a code no
// analysis phase produces on its own, so it is recognisable in a log and cannot be
// confused with a convergence failure.
#define OPS_AGREEMENT_LOST (-1000)

// THIS PROCESS'S INDEX among those that agree, and 0 when the run is not coupled.
// So `OPS_agreementRank() == 0` is the right way to ask "am I the one that should
// print the shared verdict?" in every build and every kind of run: in a parametric
// run every process answers yes, which is correct - the analyses are unrelated and
// each one must speak for itself.
int OPS_agreementRank(void);

// AGREE ON THE WORST RESULT, so every process reaches the same verdict and either
// all of them abandon the phase or none does. Failures are negative and success is
// 0, so the worst is the minimum. Returns resultHere unchanged when the run is not
// coupled, which makes every call site correct in a serial build too.
//
// MUST be called unconditionally wherever it is called at all - never only on the
// failure path. A process that reduced only when it had failed would enter this
// Allreduce while another was inside the solve's, and two mismatched collectives
// on one communicator is undefined behaviour, not a deadlock one can debug.
//
// IT ALSO CHECKS THAT THE PROCESSES ARE STILL IN STEP, at no extra cost: the
// reduce carries this process's count of agreement points alongside the result,
// once as +n and once as -n, so one MPI_INT/MPI_MIN yields both the smallest and
// the largest count. If they differ, some earlier phase returned on one process
// without agreeing, and this says so and fails the analysis instead of letting the
// run limp on with mispaired collectives. It cannot catch every case - a reduce
// that pairs with a collective of a DIFFERENT datatype is undefined behaviour and
// usually aborts before reaching here - but it does catch the quiet subset, which
// is the one that would otherwise produce wrong answers rather than no answer.
int OPS_agreeWorst(int resultHere);

// As OPS_agreeWorst, and also names the owner of the worst result: MPI_MINLOC over
// the pair (code, process) gets it out of the reduce the agreement needs anyway.
// Where several processes report the same worst code it returns the lowest-numbered
// of them, so every process names the same one and the diagnostic does not depend
// on which stderr line happens to arrive first.
//
// THE ATTRIBUTION CANNOT COME FROM THE VERDICT'S OWN REDUCE. By the time a phase
// ends, OPS_agreeWorst() has already handed every process the same negative code, so
// a MINLOC over the code alone names the lowest-numbered process instead of the one
// whose element failed. So worstRank and numFailed are reduced over whether this
// process ORIGINATED the failure rather than over the code it is now carrying, and
// numFailed goes back to meaning what it is for: one partition's element or material
// against a global problem, typically the solver.
//
// Those two extra reduces are taken on the failure path only, where the analysis is
// about to stop or sub-step anyway, and they are guarded on the ALREADY REDUCED
// verdict - the only value it is safe to guard a collective on.
//
// worstRank, numFailed and numProcesses are left untouched when the run is not
// coupled.
int OPS_agreeWorstWho(int resultHere, int &worstRank, int &numFailed, int &numProcesses);

// AGREE THAT SOMETHING HAPPENED ANYWHERE - the logical OR over the processes, and
// the flag itself when the run is not coupled. Used for domain-change detection:
// the Domain stamp is local to a process but domainChanged() is collective, so
// either all the processes call it or none does.
bool OPS_agreeAny(bool flagHere);

#endif
