# Performance on the test problems

*Last updated: 12 November 2026.*

Measured runs of SDPB on the test problems `sdpb-test-01` to `sdpb-test-08` (see [Test problems](../Readme.md#test-problems)),
one table per problem and one row per run, with comments under each table. Each problem's README gives its recommended
launch and settings, and its `reference/` folder holds the reference solution or reference iterations that a run is
checked against with `check.py`.

## Summary

Time per iteration, for the baseline (the SDPB build and settings of the Simmons-Duffin group's earlier bootstrap computations; see Setups) and for
the fastest measured run on CPUs and on GPUs, with the speedup over the baseline. The fastest run of each problem is in bold.

| Problem | N; SDP size | Baseline | Fastest on CPUs | Fastest on GPUs |
|---|---|---|---|---|
| [01](#sdpb-test-01) | 54; 0.01 GB | 0.235 s (1 node) | 0.161 s (1 node), 1.5× | **0.034 s** (1 H200), **6.9×** |
| [02](#sdpb-test-02) | 274; 0.10 GB | 1.69 s (1 node) | 1.18 s (1 node), 1.4× | **0.16 s** (1 H200), **10×** |
| [03](#sdpb-test-03) | 1,834; 2.87 GB | 43.4 s (4 nodes) | 18.8 s (4 nodes), 2.3× | **2.17 s** (1 H200), **20×** |
| [04](#sdpb-test-04) | 3,406; 9.04 GB | 70.5 s (5 nodes) | 25.8 s (5 nodes), 2.7× | **3.88 s** (1 H200), **18×** |
| [05](#sdpb-test-05) | 5,458; 8.82 GB | 95.0 s (10 nodes) | 30.3 s (10 nodes), 3.1× | **4.22 s** (2 H200s), **23×** |
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
- **Grouped Gram product on CPUs** (`--groupedSyrk`, `cpu-grouped-syrk`, #329, on top of #314). SDPB computes the
  product Q = PᵀP (P = L⁻¹B, the solve of the Schur complement blocks) exactly, in residues modulo many primes with
  double-precision BLAS (its BigInt syrk, which #314 speeds up). P is exactly zero in whole groups of columns, and in
  the same columns for many rows: on the SDPs measured (nmax = 14 to 22), the rows of P have only 6 to 19
  distinct sets of nonzero columns. With the option, the rows
  that share a set are multiplied as one group, over those columns only, so P's zeros cost no multiply-adds; Q is the
  same bit for bit, and so are the iterations. The gain follows the share of the dense multiply-adds that is left:
  0.22 on sdpb-test-05 (1.26× per iteration on 4 nodes, against the same build without the option) and 0.75 on
  sdpb-test-04 (no gain). On an nmax = 26 SDP of the group's that is not among the test problems (N = 11,002, 896
  bits), it takes 16 nodes from 105.4 to 84.0 s/iter (1.25×). The cost: no memory is saved (P is still stored dense),
  and the grouped BLAS calls, with shorter inner dimensions, run less efficiently than one dense call, which is why
  sdpb-test-04 gains nothing.
- **GPU backend** (#319–#322). It runs whole SDPB iterations on one GPU. The multiprecision matrix products are done
  exactly in residues modulo many primes on the GPU's INT8 tensor cores and recombined by the Chinese remainder theorem,
  and the blocks' Cholesky factorizations are batched by size, following Antoine Tilloy's SDPviBe. The problem must fit
  in the GPU's memory (141 GB on an H200). With #321, Q is distributed over the GPUs of a node, so a problem can use
  the memory of several (`mpirun -n G sdpb --gpu` for G GPUs); #324 extends this across nodes.
- **Faster multi-GPU Q** (`gpu-fast-q`, the seventh GPU pull request, on top of #324). Three changes to how the GPU
  backend forms Q and solves with it, all exact: the iterations are bitwise equal to #324's. (1) When Q is distributed,
  the residues of P (the solve L⁻¹B of the Schur complement blocks) are formed once per iteration and every tile of Q
  is cut from one packed result; before, each 512 × 512 tile recomputed them. (2) The matrix-vector products of the
  solves with Q are split over more threads (down to 16 terms a thread instead of 256). (3) The product Q = PᵀP skips
  P's zeros: rows of P with the same nonzero columns are multiplied as one group. These matter most on several GPUs and
  at larger N: on an nmax = 22 SDP of the group's that is not among the test problems (N = 7,990), they take two H200s
  from 13.99 to 7.65 s/iter (1.83×), and earlier runs gave 1.65× on four and 1.46× on eight. On one GPU the gain is
  about 2% (sdpb-test-03). The cost is device memory while Q is formed, up to 15 GB more per GPU on that SDP (72.3 →
  87.4 GB) and none on sdpb-test-05. `SDPB_GPU_STRUCTURED_Q=0 SDPB_GPU_Q_SINGLE_PASS=0 SDPB_GPU_GEMV_MIN_TERMS=256
  SDPB_GPU_Q_BASE=256` restores #324's code path, for example where a GPU's memory is at its limit.
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

**Fastest: one H200 with the GPU backend and the faster Q, 2.17 s/iter, 20× the baseline on 4 Expanse nodes.**
Fastest on CPUs: all CPU options, 18.8 s/iter, 2.3×.

| Setup | s/iter | Speedup | Iterations, total | Start-up | Peak memory | Check | Code | Settings |
|---|---|---|---|---|---|---|---|---|
| Baseline, Expanse 4 × 128 | 43.4 | 1× | timing: 10 iterations | 202 s (200 s) | 58.7 GB | yes: 1e-240 | 3.0.0-194-g5e92aa2f | 1024 bits |
| `--multimodularTrsm`, Expanse 4 × 128 | 22.3 | 1.9× | 269, 6,062 s (note 2) | — | — | reference | 5db9998 | 1024 bits, `--multimodularTrsm` |
| All CPU options, Expanse 4 × 128 | 18.8 | 2.3× | timing: 10 iterations | 359 s (358 s) | 64.2 GB | yes: 4e-240 | 5e725bd6 | 1024 bits, all CPU options |
| GPU backend, 1 H200 | 2.44 | 18× | 269, 658 s | 10 s | 56.7 GB (device) | yes: 1e-41 | 015f092d | 1024 bits |
| GPU backend, faster Q, 1 H200 | **2.17** | **20×** | timing: 6 iterations | 9 s | 58.0 GB (device) | yes: 5e-14 | `gpu-fast-q` 054f9bd1 | 1024 bits, `SDPB_GPU_CHUNK_BLOCKS=0` |

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
- **GPU backend, faster Q** (Setups, `gpu-fast-q`). Its 6 iterations are bitwise equal (every printed value, and y)
  to those of the branch it is built on (#324, 2f2dff3c), run on the same GPU in the same session at 2.22 s/iter; so
  on one GPU the faster Q gains about 2%, and most of the difference to the row above (2.44 s, an older build,
  another day) is not from it. The 5e-14 is the largest relative difference of the printed values from the
  reference's first 6 iterations, the same for 2f2dff3c. *Simmons-Duffin group, 2026-11-12.*

## sdpb-test-04

300 blocks (1×1 to 5×5), N = 3,406, P = 70,056, 832 bits; SDP 9.04 GB.

**Fastest: one H200 with the GPU backend, 3.88 s/iter, 18× the baseline on 5 Expanse nodes.** Fastest on
CPUs: all CPU options, 25.8 s/iter, 2.7×.

| Setup | s/iter | Speedup | Iterations, total | Start-up | Peak memory | Check | Code | Settings |
|---|---|---|---|---|---|---|---|---|
| Baseline, Expanse 5 × 128 | 70.5 | 1× | timing: 8 iterations | 212 s (208 s) | 106.9 GB | yes: 7e-162 | 3.0.0-194-g5e92aa2f | 832 bits |
| mm options, Expanse 5 × 124 | 27.8 | 2.5× | 356, 10,099 s (note 2) | — | — | reference | 2eed3fef | 832 bits, mm options |
| All CPU options, Expanse 5 × 128 | 25.8 | 2.7× | timing: 8 iterations | 7 s (none) | 109.8 GB | yes: 3e-163 | 5e725bd6 | 832 bits, all CPU options |
| All CPU options, Expanse 2 × 128 | 48.96 | 1.4× | timing: 5 iterations | 4 s (none) | — | yes: 5e-163 | `cpu-grouped-syrk` 6984360e | 832 bits, all CPU options |
| All CPU options + `--groupedSyrk`, Expanse 2 × 128 | 49.03 | 1.4× | timing: 5 iterations (2 runs) | 4–5 s (none) | — | yes: 5e-163 | `cpu-grouped-syrk` 6984360e | 832 bits, all CPU options, `--groupedSyrk` |
| GPU backend, 1 H200 | **3.88** | **18×** | 356, 1,442 s | 23 s | 127.3 GB (device) | yes: 8e-45 | 015f092d | 832 bits |

- **Baseline.** From a cold start. *Simmons-Duffin group, 2026-11-11.*
- **mm options.** The reference solution, solved to optimality in three segments, each resumed from a checkpoint. The
  Cholesky factors and solves of the Schur complement blocks and of Q by multimodular BLAS, and Q's solves on 32 ranks
  (Setups). On 124 of each node's 128 cores, where the baseline used all 128. *Simmons-Duffin group, 2026-11-03.*
- **All CPU options** (Setups), with block timings from an earlier run of the same setup, so no block-timing pass. The
  mm options, plus Q's solves and syrk reorganized for many nodes. (An earlier run stopped after 2 iterations when the
  operating system killed one MPI rank.) *Simmons-Duffin group, 2026-11-11.*
- **All CPU options, 2 nodes**, the comparison for the next row: the build of `cpu-grouped-syrk` (#329) with
  `--groupedSyrk` off, on 2 nodes, in one job with the next row's runs and with the package's
  `reference/block_timings/` for all of them. In the same job its base, `cpu-q-syrk-paths` (#314, 412644c0), ran at 48.87 s/iter, with the same printed values
  bit for bit. These runs use the stopping settings of the group's scans (error thresholds 1e-200, feasible-jump detection), which differ from `reference/sdpb-args.txt` but don't change these 5 iterations. The speedup is against the baseline's 5 nodes. *Simmons-Duffin group, 2026-11-12.*
- **All CPU options + `--groupedSyrk`** (Setups, #329), the same job and nodes as the row above: 49.18 and 48.87 s/iter
  in two runs (the row gives their mean), so no gain. Here 75% of the dense multiply-adds of Q = PᵀP remain after P's zeros are skipped, and the
  grouped BLAS calls are less efficient than one dense call. The solution x, y, the checkpoint data and every printed value are the same bit for bit
  as the row above's (only timings and the recorded options differ). *Simmons-Duffin group, 2026-11-12.*
- **GPU backend** (Setups): whole iterations on the GPU, with exact multiprecision products on its tensor cores. Solved
  in one run, with the same 356 iterations as the reference. It used 127.3 GB of the GPU's 141 GB, close to the limit
  for one H200. *Simmons-Duffin group, 2026-11-11.*

## sdpb-test-05

301 blocks (1×1 to 5×5), N = 5,458, P = 78,953, 832 bits; SDP 8.82 GB.

**Fastest: two H200s with the GPU backend and the faster Q, 4.22 s/iter, 23× the baseline on 10 Expanse nodes.**
Fastest on CPUs: the mm options, 30.3 s/iter, 3.1×.

| Setup | s/iter | Speedup | Iterations, total | Start-up | Peak memory | Check | Code | Settings |
|---|---|---|---|---|---|---|---|---|
| Baseline, Expanse 10 × 128 | 95.0 | 1× | timing: 8 iterations | — (none) | — | yes: 1e-141 | 3.0.0-194-g5e92aa2f | 832 bits |
| mm options, Expanse 10 × 128 | 30.3 | 3.1× | 384, 11,762 s (note 2) | — | — | reference | 2eed3fef | 832 bits, mm options |
| All CPU options, Expanse 10 × 128 | 33.9 | 2.8× | timing: 8 iterations | — (none) | — | yes: 1e-141 | 5e725bd6 | 832 bits, all CPU options |
| All CPU options, Expanse 4 × 128 | 56.94 | 1.7× | timing: 5 iterations (2 runs) | 7–21 s (none) | — | yes: 1e-141 | `cpu-grouped-syrk` 6984360e | 832 bits, all CPU options |
| All CPU options + `--groupedSyrk`, Expanse 4 × 128 | 45.27 | 2.1× | timing: 5 iterations (2 runs) | 7 s (none) | — | yes: 1e-141 | `cpu-grouped-syrk` 6984360e | 832 bits, all CPU options, `--groupedSyrk` |
| GPU backend, 2 H200s | 6.57 | 14× | timing: 5 iterations | 39 s | 66.3 GB (device) | yes: 3e-15 | `gpu-multinode` 2f2dff3c | 832 bits, `SDPB_GPU_SPARSE_B=1 SDPB_GPU_Q_WIDE_RESIDUAL=1` |
| GPU backend, faster Q, 2 H200s | **4.22** | **23×** | timing: 5 iterations | 38 s | 66.3 GB (device) | yes: 3e-15 | `gpu-fast-q` 054f9bd1 | 832 bits, `SDPB_GPU_SPARSE_B=1 SDPB_GPU_Q_WIDE_RESIDUAL=1` |

- **Baseline.** With the block timings of an earlier run, so no block-timing pass. From a cold
  start, that pass takes 1,676 s on these 10 nodes, and the earlier run reached no iteration within its 55-minute limit.
  The package's `reference/block_timings/` holds that pass's output, to skip it. *Simmons-Duffin group, 2026-11-11.*
- **mm options.** The reference solution, solved to optimality in three segments, each resumed from a checkpoint. The
  Cholesky factors and solves of the Schur complement blocks and of Q by multimodular BLAS, and Q's solves on 32 ranks
  (Setups). *Simmons-Duffin group, 2026-11-03.*
- **All CPU options** (Setups), in the same job as the baseline, with block timings from an earlier run. Slower than
  the mm options row on another day's nodes (note 3). *Simmons-Duffin group, 2026-11-11.*
- **All CPU options, 4 nodes**, the comparison for the next row: the build of `cpu-grouped-syrk` (#329) with
  `--groupedSyrk` off, on 4 nodes, in one job with the next row's runs (in the order off, on, base, on, off) and with
  the package's `reference/block_timings/` for all of them. Its s/iter is the mean of its two runs' (57.08 and 56.80). In the same job its
  base, `cpu-q-syrk-paths` (#314, 412644c0), ran at 56.76 s/iter, with the same printed values bit for bit. These runs use the stopping settings of the group's scans (error thresholds 1e-200, maxComplementarity 1e130, feasible-jump detection), which differ from `reference/sdpb-args.txt` but don't change these 5 iterations. The speedup is
  against the baseline's 10 nodes. *Simmons-Duffin group, 2026-11-12.*
- **All CPU options + `--groupedSyrk`** (Setups, #329), the same job and nodes as the row above: 45.41 and 45.12 s/iter
  in two runs (the row gives their mean), **1.26× the row above**. Only 22% of the dense multiply-adds of Q = PᵀP remain once the rows of P with
  the same nonzero columns are multiplied as groups over those columns. The solution x, y, the checkpoint data and every printed value are the same bit for bit as the row above's (only timings and the recorded options differ). A run in the same job with
  another valid Cholesky panel for Q (`--multimodularCholeskyQ 32`) differed in almost every file, so the comparison
  could see a change. *Simmons-Duffin group, 2026-11-12.*
- **GPU backend, 2 H200s** (Setups; #324, the branch the faster Q is built on): one node, `mpirun -n 2 --bind-to none
  sdpb --gpu`, with Q distributed over the two GPUs, B stored by its nonzeros, and the experimental residual of the
  distributed Q (`src/gpu/README.md`). The peak is the larger of the two GPUs'. *Simmons-Duffin group, 2026-11-12.*
- **GPU backend, faster Q, 2 H200s** (Setups, `gpu-fast-q`), in the same session, on the same GPUs and with the same
  launch as the row above. Its 5 iterations are bitwise equal to that row's (every printed value, and y), and 1.56×
  faster: here Q is distributed and P is 40% nonzero, so all three changes apply. The 3e-15 is the largest
  relative difference of the printed values from the reference's first 5 iterations, the same for both rows. An earlier
  whole solve with the same changes (a development build, before they were split into the pull request) took the
  reference's 384 iterations in 1,720 s, every iteration bitwise equal to the same build with the changes off.
  *Simmons-Duffin group, 2026-11-12.*

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
