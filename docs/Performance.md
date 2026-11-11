# Performance on the test problems

*Last updated: 11 November 2026.*

Measured runs of SDPB on the test problems `sdpb-test-01` to `sdpb-test-08` (see [Test problems](../Readme.md#test-problems)),
one table per problem and one row per run. Each problem's README gives its recommended launch and settings, and its
`reference/` folder holds the reference solution or reference iterations that a run is checked against with `check.py`.

**Columns.** *Code*: the commit or pull request the binary was built from. *Hardware and launch*: nodes × MPI ranks
per node, or GPUs, and the machine. *Precision and settings*: `--precision`, and any setting that differs from the
problem's `reference/sdpb-args.txt`. *Start-up*: from launch to the first iteration, including the block-timing pass
(in parentheses) that SDPB runs on CPUs when no `block_timings` file is given. *Time per iteration*: the median of
iterations 2 onward. *Iterations and total*: for a complete solve, the number of iterations and the total time; for
a timing run, "timing only: N iterations". *Peak memory*: the largest resident memory on any node (CPU runs) or on
the device (GPU runs). *Matches reference*: the outcome of `check.py`, or, for a timing run, the largest relative
difference between its iterations' printed values (objectives, gap, errors) and the reference's. *Comments*: what
was done in that run, in a sentence; longer notes are numbered under the table.

**Machines.** *Expanse* (SDSC): AMD EPYC 7742, 128 cores and 256 GB per node. *H200*: one NVIDIA H200 (141 GB) at
Caltech's Resnick cluster. *Cloud VM*: an Intel Xeon virtual machine, used only for the smallest problems.

**Code versions.** *3.0.0-194-g5e92aa2f*: the SDPB build used for production computations on Expanse, run with no
experimental options. *c5bd57e*: unmodified SDPB. Several rows come from development builds made before their code
was split into the pull requests listed in [Experimental.md](Experimental.md): 5db9998 (the code of
`cpu-mm-trsm`, #311), 2eed3fef and 5e725bd6 (the CPU options of #303, #308 and #311–#314), and 015f092d (the GPU
backend of #319–#322). A row from a published branch should name the branch and PR.

A dash means not measured. **Provisional** rows are to be repeated.

## sdpb-test-01

15 blocks (all 1×1), N = 54, P = 755, 448 bits; SDP 0.01 GB.

| Code | Hardware and launch | Precision and settings | Start-up | Time per iteration | Iterations and total | Peak memory | Matches reference | Comments | Run by | Date |
|---|---|---|---|---|---|---|---|---|---|---|
| c5bd57e | Cloud VM, 4 ranks (`mpirun`) | 448 bits | 3.0 s | 1.31 s | 290, 382 s | — | is the reference | The reference solution: solved to optimality with the settings in `sdpb-args.txt`. | Simmons-Duffin group | 2026-11-11 |
| 015f092d (GPU, #319–#322) | 1 H200, 1 rank, `--gpu` | 448 bits | < 1 s | 0.034 s | 290, 11 s | — | yes (`check.py`) | GPU backend on one GPU; the same 290 iterations as the reference. | Simmons-Duffin group | 2026-11-11 |

## sdpb-test-02

45 blocks (29 of size 1×1, 16 of size 2×2), N = 274, P = 3,758, 768 bits; SDP 0.10 GB.

| Code | Hardware and launch | Precision and settings | Start-up | Time per iteration | Iterations and total | Peak memory | Matches reference | Comments | Run by | Date |
|---|---|---|---|---|---|---|---|---|---|---|
| c5bd57e | Cloud VM, 6 ranks (`mpirun`) | 768 bits | 24 s | 13.4 s | 346, 4,651 s | — | is the reference | The reference solution: solved to optimality with the settings in `sdpb-args.txt`. | Simmons-Duffin group | 2026-11-11 |
| 015f092d (GPU, #319–#322) | 1 H200, 1 rank, `--gpu` | 768 bits | 1 s | 0.16 s | 346, 56 s | — | yes: objectives to 6e-44 | GPU backend on one GPU; the same 346 iterations as the reference. | Simmons-Duffin group | 2026-11-11 |

## sdpb-test-03

152 blocks (1×1 to 5×5), N = 1,834, P = 34,551, 1024 bits; SDP 2.87 GB.

| Code | Hardware and launch | Precision and settings | Start-up | Time per iteration | Iterations and total | Peak memory | Matches reference | Comments | Run by | Date |
|---|---|---|---|---|---|---|---|---|---|---|
| 5db9998 (#311) | Expanse, 4 × 128 (`mpirun --bind-to core`) | 1024 bits, `--multimodularTrsm` | — | 22.3 s | 269, 6,062 s (note 1) | — | is the reference | The reference solution, solved to optimality in two segments, the second resumed from a checkpoint. | Simmons-Duffin group | 2026-10-24 |
| 015f092d (GPU, #319–#322) | 1 H200, 1 rank, `--gpu` | 1024 bits | 10 s | 2.44 s | 269, 658 s | 56.7 GB (device) | yes: objectives to 1e-41 | GPU backend on one GPU, solved in one run; the same 269 iterations as the reference. | Simmons-Duffin group | 2026-11-11 |
| 3.0.0-194-g5e92aa2f | Expanse, 4 × 128 (`srun`) | 1024 bits | 202 s (200 s) | 43.4 s | timing only: 10 iterations | 58.7 GB | yes: 1e-240 | Production build with no experimental options, from a cold start (block-timing pass included). | Simmons-Duffin group | 2026-11-11 |
| 5e725bd6 (#303, #308, #311–#314) | Expanse, 4 × 128 (`srun`) | 1024 bits, note 2 | 359 s (358 s) | 18.8 s | timing only: 10 iterations | 64.2 GB | yes: 4e-240 | The CPU options of the listed PRs on (note 2), same nodes and launch as the row above; 2.3× faster per iteration. | Simmons-Duffin group | 2026-11-11 |

## sdpb-test-04

300 blocks (1×1 to 5×5), N = 3,406, P = 70,056, 832 bits; SDP 9.04 GB.

| Code | Hardware and launch | Precision and settings | Start-up | Time per iteration | Iterations and total | Peak memory | Matches reference | Comments | Run by | Date |
|---|---|---|---|---|---|---|---|---|---|---|
| 2eed3fef (#311–#313) | Expanse, 5 × 124 (`srun`) | 832 bits, `--multimodularTrsm --multimodularCholeskyQ 64 --qSolveRanks 32` | — | 27.8 s | 356, 10,099 s (note 1) | — | is the reference | The reference solution, solved to optimality in three segments, each resumed from a checkpoint. | Simmons-Duffin group | 2026-11-03 |
| 015f092d (GPU, #319–#322) | 1 H200, 1 rank, `--gpu` | 832 bits | 23 s | 3.88 s | 356, 1,442 s | 127.3 GB (device) | yes: objectives to 8e-45 | GPU backend on one GPU, solved in one run; the same 356 iterations as the reference. | Simmons-Duffin group | 2026-11-11 |
| 3.0.0-194-g5e92aa2f | Expanse, 5 × 128 (`srun`) | 832 bits | 212 s (208 s) | 70.5 s | timing only: 8 iterations | 106.9 GB | yes: 7e-162 | Production build with no experimental options, from a cold start. | Simmons-Duffin group | 2026-11-11 |
| 5e725bd6 (#303, #308, #311–#314) **provisional** | Expanse, 5 × 128 (`srun`) | 832 bits, note 2 | 8 s (none) | 29.7 s (1 iteration) | stopped after 2 iterations (note 3) | 103.5 GB | yes: 1e-163 (2 iterations) | The CPU options of the listed PRs on, with block timings reused from an earlier run (no block-timing pass); one MPI rank was killed in iteration 3. To be repeated. | Simmons-Duffin group | 2026-11-11 |

## sdpb-test-05

301 blocks (1×1 to 5×5), N = 5,458, P = 78,953, 832 bits; SDP 8.82 GB.

| Code | Hardware and launch | Precision and settings | Start-up | Time per iteration | Iterations and total | Peak memory | Matches reference | Comments | Run by | Date |
|---|---|---|---|---|---|---|---|---|---|---|
| 2eed3fef (#311–#313) | Expanse, 10 × 128 (`srun`) | 832 bits, `--multimodularTrsm --multimodularCholeskyQ 64 --qSolveRanks 32` | — | 30.3 s | 384, 11,762 s (note 1) | — | is the reference | The reference solution, solved to optimality in three segments, each resumed from a checkpoint. | Simmons-Duffin group | 2026-11-03 |
| 3.0.0-194-g5e92aa2f | Expanse, 10 × 128 (`srun`) | 832 bits | 1,680 s (1,676 s) | — | none within the 55-minute limit | 114.7 GB | — | Production build from a cold start: its block-timing pass used most of the time limit. The package's `reference/block_timings/` holds this pass's output, to skip it in later runs. | Simmons-Duffin group | 2026-11-11 |
| 5e725bd6 (#303, #308, #311–#314) **provisional** | Expanse, 10 × 128 (`srun`) | 832 bits, note 2 | 352 s (350 s) | 74.8 s | timing only: 3 iterations | 118.8 GB | yes: 1e-141 | The CPU options of the listed PRs on, from a cold start. 2.5× slower per iteration than the reference row with similar options on as many nodes (note 4). To be repeated. | Simmons-Duffin group | 2026-11-11 |

## sdpb-test-06, sdpb-test-07

To be added.

## sdpb-test-08

384 blocks (1×1 to 8×8), N = 3,047, P = 84,576, 768 bits; SDP 3.28 GB. A timing problem: its reference is the first
8 iterations of the production build (first row), not a complete solve.

| Code | Hardware and launch | Precision and settings | Start-up | Time per iteration | Iterations and total | Peak memory | Matches reference | Comments | Run by | Date |
|---|---|---|---|---|---|---|---|---|---|---|
| 3.0.0-194-g5e92aa2f | Expanse, 2 × 128 (`srun`) | 768 bits | 378 s (378 s) | 109.0 s | timing only: 8 iterations | 208.8 GB | is the reference | Production build with no experimental options, from a cold start; its 8 iterations are the reference iterations. | Simmons-Duffin group | 2026-11-11 |
| 5e725bd6 (#303, #308, #311–#314) | Expanse, 2 × 128 (`srun`) | 768 bits, note 2 | 161 s (160 s) | 41.4 s | timing only: 8 iterations | — | yes: 1.2e-189 | The CPU options of the listed PRs on, from a cold start, in the same job as the row above; 2.6× faster per iteration. | Simmons-Duffin group | 2026-11-11 |

## SDPA-format problems (`sdpb-test-sdpa-NN`)

To be added: `.dat-s` problems of increasing size, run with `sdpb`'s SDPA mode
(see [Solving SDPA problems](Usage.md#solving-sdpa-problems)), on CPU and GPU.

## Notes

1. The total is the sum of the iterations' times over the segments; the start-up of each segment is not included.
2. `--multimodularTrsm --multimodularCholeskyQ 64 --qSolveRanks 32`, with the node-shared diagonal blocks and
   factor, the solve cache and the faster syrk paths on. In this development build the last four were set by
   environment variables; on the published branches they are the options `--qCholeskyNodeDiag`,
   `--qCholeskyNodeFactor`, `--qSolveCache` (#313) and `--batchCrt --batchResidues --residuesNoClear
   --residuesStarVR` (#314).
3. One MPI rank was killed by the operating system during iteration 3 on one node, and the job stopped. The two
   iterations it completed match the reference.
4. Not understood yet. The two runs used different nodes; a node or filesystem effect has not been ruled out, and
   the run will be repeated on other nodes.

## How to add a result

1. Run the problem as its README describes, on your hardware and with your build.
2. Run `python3 check.py <your outDir>` (or the README's command for timing problems) and note what it reports.
3. Open a pull request that adds one row to the problem's table, with every column above filled in: the commit or
   PR the binary was built from, the hardware and launch, the precision and any settings that differ from
   `sdpb-args.txt`, the start-up time, the time per iteration, the iterations and total time (or "timing only: N
   iterations"), the peak memory, the result of `check.py`, a **comment** saying what was done (which options or
   PRs were on, what changed from the row it compares with, anything unusual such as reused block timings or an
   excluded node, and the outcome if it isn't a plain match), your name and the date. Longer notes go under the
   table, numbered.
4. One row per configuration: a run with different options, launch or hardware gets its own row.
