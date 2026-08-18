# system ClusterPardiso (OpenSeesMP)

Distributed direct solver on MKL `cluster_sparse_solver` (Cluster PARDISO).
Each MPI rank owns a contiguous block of equations and hands its rows of the
assembled matrix to the solver in CSR form (`iparm[39]=2` distributed input);
element contributions to rows owned by other ranks are exchanged once per
`setSize` (pattern) and at each factorization (values).

## Usage

```tcl
# OpenSeesMP only (mpiexec -n P OpenSeesMP.exe script.tcl)
constraints Transformation
numberer   ParallelPlain        ;# REQUIRED for np > 1 (equation ownership)
system     ClusterPardiso       ;# defaults: safe for classic FEM analyses
```

The defaults are the MKL-recommended robust settings for real unsymmetric
matrices (mtype 11): METIS ordering, scaling + weighted matching on, double
precision, in-core. They are the configuration validated to 10 digits against
`system Mumps` on linear static, eigen and nonlinear (J2) benchmarks. Users
who do not pass any flag get exactly that.

## Flags

Indices are the **0-based C indices of the MKL Developer Reference**
(`iparm[i]` as written in C code). Beware: some MKL manual pages use 1-based
Fortran numbering — `iparm(28)` in Fortran is `iparm[27]` in C.

| Flag | Maps to | Default | Effect |
|---|---|---|---|
| `-spd` | `mtype=2` | 11 | SPD matrix, Cholesky LL^T on the upper triangle only: ~half memory and factorization flops, backsolves read half the factor. Measured on a 620k-eqn wall: static 1.7x, eigen 1.5x, nonlinear 1.3x faster than the mtype-11 default (and faster than Mumps on all three). Fails with error -4 if the matrix is not positive definite |
| `-sym` | `mtype=-2` | 11 | symmetric indefinite, LDL^T with Bunch-Kaufman pivoting: same storage benefits as `-spd`, tolerates indefinite tangents (limit points, softening) |
| `-ordering <0\|2\|3>` | `iparm[1]` | 2 | fill-reducing ordering: 0 = minimum degree, 2 = METIS nested dissection, 3 = OpenMP-parallel METIS (useful with threaded MKL). Advisory: with np > 1 the cluster solver appears to ignore 0 (measured: identical result digits to METIS) |
| `-refine <n>` | `iparm[7]` | 0 | up to n steps of iterative refinement per solve (measured: `-refine 2` halved the residual on the validation wall) |
| `-noscale` | `iparm[10]=0`, `iparm[12]=0` | on | disable scaling + weighted matching; safe and slightly faster for SPD-like FE matrices. Implied by `-spd`/`-sym` (unsymmetric-only features) |
| `-msglvl` | `msglvl=1` | 0 | print MKL statistics (reordering, memory, flops) |
| `-iparm <i> <v>` | `iparm[i]=v` | — | raw override, repeatable, for experimentation without recompiling |

With `-spd`/`-sym` the assembly keeps the full pattern (its cost is
negligible); only the upper triangle is handed to MKL and exchanged across
ranks. Structural check: csrNnz summed over ranks = (full nnz + n) / 2.

Example:

```tcl
system ClusterPardiso -noscale -refine 1
system ClusterPardiso -iparm 23 1        ;# raw: two-level factorization scheduling
```

Note on out-of-core: `iparm[59]=2` was tried and REFUSED by
`cluster_sparse_solver` (error -1, inconsistent input, at phase 11 with
np = 2) — OOC appears to be a serial-PARDISO feature not supported by the
cluster version, so there is deliberately no `-ooc` flag. If you want to
re-test on a newer MKL: `system ClusterPardiso -iparm 59 2`.

## Reserved iparm indices

`-iparm` refuses these with an error — they are the contract between the
OpenSees integration and MKL, and overriding them breaks the solver in
non-obvious ways:

| Index | Owned meaning |
|---|---|
| 0 | "user iparm" switch (always 1) |
| 34 | 0/1-based indexing (integration uses 1-based ia/ja) |
| 39 | distributed input format (=2, row-block CSR) |
| 40, 41 | this rank's owned row range (set from EquationPartition) |

Flags must appear on the `system` line: options are refused after the first
analysis has started.

## Debugging

`OPS_CPARDISO_DEBUG=1` (environment variable) enables:
- MKL input checking (`iparm[26]=1`),
- a residual print `||A*x-b|| / ||b||` after every solve, computed
  independently from both the assembled COO and the CSR handed to MKL
  (layer separation: assembly vs mapping vs solver).

Solver errors are printed with the MKL error code and a short explanation
(e.g. `-2: not enough memory`, `-4: zero pivot ... singular matrix?`).

## Known constraints and planned flags

- `numberer ParallelPlain` is mandatory for np > 1; other numberers give a
  clean error (no equation-ownership blocks).
- Phase 11 (analysis) runs at the FIRST solve, not at `setSize`: with
  scaling/matching on, MKL's analysis reads the numerical values, and at
  `setSize` time the matrix is still zero. Do not "optimize" this back.
- The `OPS_CPARDISO_DEBUG` CSR residual check reports `n/a(sym)` with
  `-spd`/`-sym`: the upper-triangle CSR alone cannot reproduce full row
  sums. The COO residual stays exact (assembly keeps the full pattern).
- Planned, deliberately not exposed yet:
  - `-single` — FP32 factorization (`iparm[27]=1`); MKL then requires
    float a/b/x buffers, so it needs conversion buffers in the solver.
  - `-krylov` — `iparm[3]` CGS/CG preconditioned by a stale factorization;
    experimental, for nonlinear analyses.
