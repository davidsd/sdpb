# Performance on the test problems

*Last updated: 11 November 2026.*

Measured runs of SDPB on the test problems `sdpb-test-01` to `sdpb-test-08` (see [Test problems](../Readme.md#test-problems)),
one table per problem and one row per run, with comments under each table. Each problem's README gives its recommended
launch and settings, and its `reference/` folder holds the reference solution or reference iterations that a run is
checked against with `check.py`.

## Summary

Time per iteration, for the baseline (the SDPB build and settings of the Simmons-Duffin group's earlier bootstrap computations; see Setups) and for
the fastest measured run on CPUs and on a GPU, with the speedup over the baseline. The fastest run of each problem is in bold.

| Problem | N; SDP size | Baseline | Fastest on CPUs | Fastest on a GPU |
|---|---|---|---|---|
| [01](#sdpb-test-01) | 54; 0.01 GB | 0.235 s (1 node) | 0.161 s (1 node), 1.5× | **0.034 s** (1 H200), **6.9×** |
| [02](#sdpb-test-02) | 274; 0.10 GB | 1.69 s (1 node) | 1.18 s (1 node), 1.4× | **0.16 s** (1 H200), **10×** |
| [03](#sdpb-test-03) | 1,834; 2.87 GB | 43.4 s (4 nodes) | 18.8 s (4 nodes), 2.3× | **2.44 s** (1 H200), **18×** |
| [04](#sdpb-test-04) | 3,406; 9.04 GB | 70.5 s (5 nodes) | 25.8 s (5 nodes), 2.7× | **3.88 s** (1 H200), **18×** |
| [05](#sdpb-test-05) | 5,458; 8.82 GB | 95.0 s (10 nodes) | **30.3 s** (10 nodes), **3.1×** | — |
| [08](#sdpb-test-08) | 3,047; 3.28 GB | 109.0 s (2 nodes) | **41.4 s** (2 nodes), **2.6×** | — |

Nodes are Expanse nodes (128 cores each). A dash means not measured; 06 and 07 are to be added.

## Setups

How each code in the tables gets its speed. The pull requests are described in [Experimental.md](Experimental.md).

- **Baseline** (3.0.0-194-g5e92aa2f, on the `spectrum-mpsolve` branch of SDPB, April 2025): the SDPB build and solver
  settings that the Simmons-Duffin group used for its bootstrap computations in 2025, with no experimental options.
  Every problem's speedups are relative to it.
- **Multimodular CPU options** (#311–#313, with #303 and #314). SDPB's most expensive steps are products and
  factorizations of matrices of multiprecision (BigFloat) numbers. `--multimodularTrsm` (#311) computes each Schur
  complement block's Cholesky factor L and the solve L⁻¹B exactly by multimodular arithmetic: fixed-point integers are
  taken from the BigFloat values, reduced modulo about 100 primes, multiplied with double-precision BLAS (`dgemm`), and
  recombined by the Chinese remainder theorem. Only small diagonal blocks stay in BigFloat. `--multimodularCholeskyQ 64`
  (#312) does the same for the Cholesky decomposition of Q, the N × N matrix of the Schur complement. `--qSolveRanks 32`
  (#313) does the solves with Q's factor on 32 ranks instead of all of them, which cuts communication on many nodes.
  These three options are called the **mm options** in the tables. **All CPU options** adds the node-shared diagonal
  blocks and factor of Q and a solve cache (#313), the faster paths of Q's BigInt syrk (#314), and a
  `restore_and_reduce` that visits only each rank's own entries of Q (#303).
  The cost: the results are not bitwise equal to `master`'s. The multimodular steps differ within a proven error bound
  far below the working precision, and `--qSolveRanks` changes the iterates at working precision. Each row's Check
  column says how closely its iterations match the reference.
- **GPU backend** (#319–#322). It runs whole SDPB iterations on one GPU. The multiprecision matrix products are done
  exactly in residues modulo many primes on the GPU's INT8 tensor cores and recombined by the Chinese remainder theorem,
  and the blocks' Cholesky factorizations are batched by size, following Antoine Tilloy's SDPviBe. The problem must fit
  in the GPU's memory (141 GB on an H200); distributing Q over several GPUs (#321) has not been timed here.
- **Development builds.** Several rows come from builds made before their code was split into the pull requests above:
  5db9998 (the code of #311), 2eed3fef (#311–#313), 5e725bd6 (#303, #308, #311–#314) and 015f092d (#319–#322). In
  5e725bd6, the options of #313 other than `--qSolveRanks` and the four of #314 were set by environment variables; on the
  published branches they are `--qCholeskyNodeDiag`, `--qCholeskyNodeFactor`, `--qSolveCache` (#313) and `--batchCrt
  --batchResidues --residuesNoClear --residuesStarVR` (#314). A row from a published branch should name the branch and PR.

## Columns and machines

*Setup*: the code (above) and the hardware: Expanse nodes × MPI ranks per node, or GPUs. *s/iter*: the time per
iteration, the median of iterations 2 onward; the problem's fastest run and its speedup are in bold. *Speedup*: the
baseline row's s/iter divided by this row's. *Iterations, total*: for a complete solve, the number of
iterations and the total time; for a timing run, "timing: N iterations". *Start-up*: from launch to the first
iteration, with the block-timing pass that SDPB runs on CPUs when no `block_timings` file is given in parentheses.
*Peak memory*: the largest resident memory on any node (CPU runs) or on the device (GPU runs). *Check*: the outcome of
`check.py`, or, for a timing run, the largest relative difference between its iterations' printed values (objectives,
gap, errors) and the reference's. *Code*: the commit the binary was built from. *Settings*: `--precision`, and the
options that differ from the problem's `reference/sdpb-args.txt`. Under each table, a comment per row says what was
done, how the speedup was obtained, who ran it and when.

*Expanse* (SDSC): AMD EPYC 7742, 128 cores and 256 GB per node; runs are launched with `srun` unless a comment says
otherwise. *H200*: one NVIDIA H200 (141 GB) at Caltech's Resnick cluster.

A dash means not measured.

## sdpb-test-01

15 blocks (all 1×1), N = 54, P = 755, 448 bits; SDP 0.01 GB. Its reference solution was computed with SDPB's `master`
(3.1.0-27-gc5bd57e9).

**Fastest: one H200 with the GPU backend, 0.034 s/iter, 6.9× the baseline on one Expanse node.**

| Setup | s/iter | Speedup | Iterations, total | Start-up | Peak memory | Check | Code | Settings |
|---|---|---|---|---|---|---|---|---|
| Baseline, Expanse 1 × 128 | 0.235 | 1× | 290, — | 8 s (8 s) | — | yes: 2e-40 | 3.0.0-194-g5e92aa2f | 448 bits |
| All CPU options, Expanse 1 × 128 | 0.161 | 1.5× | 290, — | 1 s (none) | 15.0 GB | yes: 2e-40 | 5e725bd6 | 448 bits, all CPU options |
| GPU backend, 1 H200 | **0.034** | **6.9×** | 290, 11 s | < 1 s | — | yes (`check.py`) | 015f092d | 448 bits |

- **Baseline.** Solved to optimality from a cold start, on an `sdp/` converted again on Expanse
  (note 1). The same optimum as the reference. *Simmons-Duffin group, 2026-11-11.*
- **All CPU options** (Setups), on the same node count as the baseline and in the same job. The speedup is
  smaller than on the larger problems. *Simmons-Duffin group, 2026-11-11.*
- **GPU backend** (Setups): whole iterations on the GPU, with exact multiprecision products on its tensor cores. The
  same 290 iterations as the reference. *Simmons-Duffin group, 2026-11-11.*

## sdpb-test-02

45 blocks (29 of size 1×1, 16 of size 2×2), N = 274, P = 3,758, 768 bits; SDP 0.10 GB. Its reference solution was
computed with SDPB's `master` (3.1.0-27-gc5bd57e9).

**Fastest: one H200 with the GPU backend, 0.16 s/iter, 10× the baseline on one Expanse node.**

| Setup | s/iter | Speedup | Iterations, total | Start-up | Peak memory | Check | Code | Settings |
|---|---|---|---|---|---|---|---|---|
| Baseline, Expanse 1 × 128 | 1.69 | 1× | 346, — | 12 s (12 s) | — | yes: 1e-84 | 3.0.0-194-g5e92aa2f | 768 bits |
| All CPU options, Expanse 1 × 128 | 1.18 | 1.4× | 346, — | 5 s (none) | — | yes: 3e-84 | 5e725bd6 | 768 bits, all CPU options |
| GPU backend, 1 H200 | **0.16** | **10×** | 346, 56 s | 1 s | — | yes: 6e-44 | 015f092d | 768 bits |

- **Baseline.** Solved to optimality from a cold start, on an `sdp/` converted again on Expanse
  (note 1). The same optimum as the reference. *Simmons-Duffin group, 2026-11-11.*
- **All CPU options** (Setups), on the same node count as the baseline and in the same job. *Simmons-Duffin group,
  2026-11-11.*
- **GPU backend** (Setups): whole iterations on the GPU, with exact multiprecision products on its tensor cores. The
  same 346 iterations as the reference, and the same optimum to 6e-44. *Simmons-Duffin group, 2026-11-11.*

## sdpb-test-03

152 blocks (1×1 to 5×5), N = 1,834, P = 34,551, 1024 bits; SDP 2.87 GB.

**Fastest: one H200 with the GPU backend, 2.44 s/iter, 18× the baseline on 4 Expanse nodes.** Fastest on
CPUs: all CPU options, 18.8 s/iter, 2.3×.

| Setup | s/iter | Speedup | Iterations, total | Start-up | Peak memory | Check | Code | Settings |
|---|---|---|---|---|---|---|---|---|
| Baseline, Expanse 4 × 128 | 43.4 | 1× | timing: 10 iterations | 202 s (200 s) | 58.7 GB | yes: 1e-240 | 3.0.0-194-g5e92aa2f | 1024 bits |
| `--multimodularTrsm`, Expanse 4 × 128 | 22.3 | 1.9× | 269, 6,062 s (note 2) | — | — | reference | 5db9998 | 1024 bits, `--multimodularTrsm` |
| All CPU options, Expanse 4 × 128 | 18.8 | 2.3× | timing: 10 iterations | 359 s (358 s) | 64.2 GB | yes: 4e-240 | 5e725bd6 | 1024 bits, all CPU options |
| GPU backend, 1 H200 | **2.44** | **18×** | 269, 658 s | 10 s | 56.7 GB (device) | yes: 1e-41 | 015f092d | 1024 bits |

- **Baseline.** From a cold start (block-timing pass included). *Simmons-Duffin group, 2026-11-11.*
- **`--multimodularTrsm`.** The reference solution, solved to optimality in two segments, the second resumed from a
  checkpoint. Only #311 is on: the Cholesky factor of each Schur complement block and the solve L⁻¹B are computed by
  multimodular BLAS instead of BigFloat arithmetic (Setups). Launched with `mpirun --bind-to core`.
  *Simmons-Duffin group, 2026-10-24.*
- **All CPU options** (Setups), on the same nodes and launch as the baseline and in the same job: the multimodular
  factors and solves of the row above, the same for Q, and Q's solves and syrk reorganized for many nodes. The cost
  here: 5.5 GB more peak memory per node, and a longer block-timing pass. *Simmons-Duffin group, 2026-11-11.*
- **GPU backend** (Setups): whole iterations on the GPU, with exact multiprecision products on its tensor cores. Solved
  in one run, with the same 269 iterations as the reference. It used 56.7 GB of the GPU's 141 GB.
  *Simmons-Duffin group, 2026-11-11.*

## sdpb-test-04

300 blocks (1×1 to 5×5), N = 3,406, P = 70,056, 832 bits; SDP 9.04 GB.

**Fastest: one H200 with the GPU backend, 3.88 s/iter, 18× the baseline on 5 Expanse nodes.** Fastest on
CPUs: all CPU options, 25.8 s/iter, 2.7×.

| Setup | s/iter | Speedup | Iterations, total | Start-up | Peak memory | Check | Code | Settings |
|---|---|---|---|---|---|---|---|---|
| Baseline, Expanse 5 × 128 | 70.5 | 1× | timing: 8 iterations | 212 s (208 s) | 106.9 GB | yes: 7e-162 | 3.0.0-194-g5e92aa2f | 832 bits |
| mm options, Expanse 5 × 124 | 27.8 | 2.5× | 356, 10,099 s (note 2) | — | — | reference | 2eed3fef | 832 bits, mm options |
| All CPU options, Expanse 5 × 128 | 25.8 | 2.7× | timing: 8 iterations | 7 s (none) | 109.8 GB | yes: 3e-163 | 5e725bd6 | 832 bits, all CPU options |
| GPU backend, 1 H200 | **3.88** | **18×** | 356, 1,442 s | 23 s | 127.3 GB (device) | yes: 8e-45 | 015f092d | 832 bits |

- **Baseline.** From a cold start. *Simmons-Duffin group, 2026-11-11.*
- **mm options.** The reference solution, solved to optimality in three segments, each resumed from a checkpoint. The
  Cholesky factors and solves of the Schur complement blocks and of Q by multimodular BLAS, and Q's solves on 32 ranks
  (Setups). On 124 of each node's 128 cores, where the baseline used all 128. *Simmons-Duffin group, 2026-11-03.*
- **All CPU options** (Setups), with block timings from an earlier run of the same setup, so no block-timing pass. The
  mm options, plus Q's solves and syrk reorganized for many nodes. (An earlier run stopped after 2 iterations when the
  operating system killed one MPI rank.) *Simmons-Duffin group, 2026-11-11.*
- **GPU backend** (Setups): whole iterations on the GPU, with exact multiprecision products on its tensor cores. Solved
  in one run, with the same 356 iterations as the reference. It used 127.3 GB of the GPU's 141 GB, close to the limit
  for one H200. *Simmons-Duffin group, 2026-11-11.*

## sdpb-test-05

301 blocks (1×1 to 5×5), N = 5,458, P = 78,953, 832 bits; SDP 8.82 GB.

**Fastest: the mm options on 10 Expanse nodes, 30.3 s/iter, 3.1× the baseline on the same nodes.** Not yet
run on a GPU.

| Setup | s/iter | Speedup | Iterations, total | Start-up | Peak memory | Check | Code | Settings |
|---|---|---|---|---|---|---|---|---|
| Baseline, Expanse 10 × 128 | 95.0 | 1× | timing: 8 iterations | — (none) | — | yes: 1e-141 | 3.0.0-194-g5e92aa2f | 832 bits |
| mm options, Expanse 10 × 128 | **30.3** | **3.1×** | 384, 11,762 s (note 2) | — | — | reference | 2eed3fef | 832 bits, mm options |
| All CPU options, Expanse 10 × 128 | 33.9 | 2.8× | timing: 8 iterations | — (none) | — | yes: 1e-141 | 5e725bd6 | 832 bits, all CPU options |

- **Baseline.** With the block timings of an earlier run, so no block-timing pass. From a cold
  start, that pass takes 1,676 s on these 10 nodes, and the earlier run reached no iteration within its 55-minute limit.
  The package's `reference/block_timings/` holds that pass's output, to skip it. *Simmons-Duffin group, 2026-11-11.*
- **mm options.** The reference solution, solved to optimality in three segments, each resumed from a checkpoint. The
  Cholesky factors and solves of the Schur complement blocks and of Q by multimodular BLAS, and Q's solves on 32 ranks
  (Setups). *Simmons-Duffin group, 2026-11-03.*
- **All CPU options** (Setups), in the same job as the baseline, with block timings from an earlier run. Slower than
  the mm options row on another day's nodes (note 3). *Simmons-Duffin group, 2026-11-11.*

## sdpb-test-06, sdpb-test-07

To be added.

## sdpb-test-08

384 blocks (1×1 to 8×8), N = 3,047, P = 84,576, 768 bits; SDP 3.28 GB. A timing problem: its reference is the first
8 iterations of the baseline (first row), not a complete solve.

**Fastest: all CPU options on 2 Expanse nodes, 41.4 s/iter, 2.6× the baseline on the same nodes.** Not yet
run on a GPU.

| Setup | s/iter | Speedup | Iterations, total | Start-up | Peak memory | Check | Code | Settings |
|---|---|---|---|---|---|---|---|---|
| Baseline, Expanse 2 × 128 | 109.0 | 1× | timing: 8 iterations | 378 s (378 s) | 208.8 GB | reference | 3.0.0-194-g5e92aa2f | 768 bits |
| All CPU options, Expanse 2 × 128 | **41.4** | **2.6×** | timing: 8 iterations | 161 s (160 s) | — | yes: 1.2e-189 | 5e725bd6 | 768 bits, all CPU options |

- **Baseline.** From a cold start; its 8 iterations are the reference iterations.
  *Simmons-Duffin group, 2026-11-11.*
- **All CPU options** (Setups), from a cold start, in the same job and on the same nodes as the baseline. The Cholesky
  factors and solves of the Schur complement blocks and of Q are done by multimodular BLAS instead of BigFloat
  arithmetic. Its block-timing pass also took less than half as long as the baseline's. *Simmons-Duffin group,
  2026-11-11.*

## SDPA-format problems (`sdpb-test-sdpa-NN`)

To be added: `.dat-s` problems of increasing size, run with `sdpb`'s SDPA mode
(see [Solving SDPA problems](Usage.md#solving-sdpa-problems)), on CPU and GPU.

## Notes

1. The Expanse runs of 01 and 02 read an `sdp/` converted from the package's PMP on Expanse, because the shipped
   `sdp/` was written with a newer Boost than the Expanse builds read.
2. The total is the sum of the iterations' times over the segments; the start-up of each segment is not included.
3. Not understood yet. The two runs used different nodes on different days; per-iteration times of these options have
   been seen to move by 5–10% between allocations. A first run of all CPU options, on other nodes and from a cold
   start, gave 74.8 s/iter over 3 iterations; this rerun gave 33.9.

## How to add a result

1. Run the problem as its README describes, on your hardware and with your build.
2. Run `python3 check.py <your outDir>` (or the README's command for timing problems) and note what it reports.
3. Open a pull request that adds one row to the problem's table, with every column filled in: a short setup (code and
   hardware), the time per iteration, the speedup over the table's baseline row, the iterations and total time (or
   "timing: N iterations"), the start-up time, the peak memory, the result of `check.py`, the commit or PR the binary
   was built from, and the precision and any settings that differ from `sdpb-args.txt`.
4. Add a comment under the table, starting with your row's setup in bold. Explain how the speedup was obtained, in a
   few sentences: what was changed (algorithm, kernel, data layout, communication, precision, launch or settings), why
   it is faster, and what it costs or gives up (memory, accuracy, generality), with a link to the PR or commit. Say
   also what was unusual about the run, such as reused block timings or an excluded node. End with your name and the
   date. If your run is the new fastest, update the line above the table, the bold and the summary.
5. One row per configuration: a run with different options, launch or hardware gets its own row.
