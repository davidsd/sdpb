# GPU backend for SDPB (in development)

**The GPU backend is based on [SDPviBe](https://github.com/AntoineTilloy/SDPviBe) by Antoine Tilloy, the
prototype it follows.** SDPviBe re-implements SDPB 3.1's interior-point method on GPUs; this backend brings its
algorithms into SDPB itself, in C++/CUDA, so that SDPB stays one codebase. Taken from SDPviBe (details and file
mapping in [docs/GPU_backend_credits.md](../../docs/GPU_backend_credits.md)):

- exact multiprecision matrix products in a residue number system (RNS), with the products on INT8 tensor
  cores and the reconstruction (CRT) as an INT8 product with the digits of M/m_i;
- the Gram matrix Q = PᵀP accumulated in residue form over chunks, with one reconstruction per iteration;
- batched, equilibrated fixed-point Cholesky and triangular inverses over SDP blocks of equal size;
- step lengths from float64 eigenvalues with an error bound and multiprecision refinement;
- chunked streaming of per-block arrays between host and device (planned), and NCCL collectives
  (the experimental distributed-Q backend below).

Device-float products first round operands to the chosen per-column
alignment. RNS accumulation is exact for those aligned integers; the result
is then rounded to a floating mantissa. The real-Q accuracy/cost tradeoff and
explicit options are documented in [the regression recipe](test/q_accuracy.md).

Float sizes: `rns/Limb_Sizes.hxx` lists the L with kernels compiled for them (even L from 24 to 40: --precision up to 1246 bits, rounded up to L = 24 below 704); adding a size is one line there and costs compile time. Checks take `--L`.

## What is here (products and block factorization)

- `rns/Rns_Basis.{hxx,cxx}`: primes in (2^13, 2^14), residue and CRT tables (host, GMP).
- `rns/Gram_Accumulator.{hxx,cu}`: exact Q = Σ_chunks AᵀA for big-integer A (balanced base-256 digits),
  residues and products through cuBLASLt INT8 GEMMs (wgmma on Hopper), three INT8 products per prime
  (Karatsuba; SDPviBe uses four), the lower triangle cut into recursive rectangles, uint16 residues of the
  packed lower triangle kept across chunks, CRT once in `finish()`.
- `rns/Mpf_Gram.{hxx,cu}`: the same on device floats (SDPviBe's format, column-major): column maxima,
  residues straight from the 32-bit limbs, exact Gram, each entry of Q rounded once to a float.
  `add_structured` takes the product by groups of rows with the same nonzero columns (SDPB's P = L^-1 B is zero in
  whole groups of columns): the same Q bit for bit, with the multiply-adds of the nonzero part only.
- `rns/Rns_Gemm.{hxx,cu}`: exact batched products C_b = A_bᵀB_b over equal-size blocks (L⁻¹B with an
  explicit inverse, pairings, refinement), column panels of B, one CRT per panel.
- `rns/rns_kernels.cuh`: shared kernels: residue finish, CRT, fused residues and fused CRT (WMMA INT8).
- `linalg/basecase.cuh`: SDPviBe's CUDA-core base case (batched 32×32 multiprecision Cholesky and
  triangular inverse in wide fixed point), ported to templates.
- `linalg/Block_Factor.{hxx,cu}`: recursive Cholesky and inverse, with whole-matrix
  power-of-two scaling and one panel residual correction by default. The original
  unrefined recursion is available for regression comparisons. Adaptive correction
  precision is optional. See [the real-Q accuracy checks](test/q_accuracy.md) for
  the tested scope, independent references and experimental options.
- `Solver_Backend.hxx`: the interface between SDPB's driver and the backend (design sketch).
- `test/*_check.cu`: exactness/accuracy against GMP (independent schoolbook products), and timing.
- `bench/sdpvibe_*_bench.py`: SDPviBe's own GramAcc, gemm, Chol, solve_panels on the same shapes.

- `solver/Gpu_Solver.{hxx,cu}`: SDPB's arithmetic per iteration on one GPU, below the driver `SDP_Solver::run`:
  blocks batched by dimension with inert padding, as SDPviBe's solver; residues, Schur blocks, Q, search directions and
  step lengths. `solver/Mpf_Ops.{hxx,cu}`: element kernels, assemblies, direct multiprecision matrix-vector products
  and reductions. `src/sdp_solve/SDP_Solver/run/gpu/Gpu_Compute.{hxx,cxx}`: the conversions between `El::BigFloat`
  and device floats, and what `SDP_Solver::run` calls.

## sdpb --gpu

`./waf configure --enable-gpu --nvcc="nvcc --compiler-bindir <gcc 13 bin>" --gpu-arch=90 --cuda-dir=<CUDA> \
--gpu-gmp-include=<GMP include> <the usual options>`, then `./waf --targets=sdpb`. nvcc's host compiler may
differ from SDPB's (the interface header has no CUDA or Elemental types). Without `--enable-gpu` the CPU build is unchanged.
With CMake: `cmake -B build -DSDPB_ENABLE_GPU=ON -DCMAKE_CUDA_HOST_COMPILER=<gcc 13 bin>/g++ <the usual options>`,
then `cmake --build build --target sdpb` (CUDA architectures: `-DCMAKE_CUDA_ARCHITECTURES`, `90-real` by default).
The option is off by default.

With `--controlDir` the process stays resident between calls, SDP and iterate on the device (`src/sdpb/resident/Readme.md`). `mpirun -n 1 sdpb --gpu <options>` runs the iterations on one GPU (one MPI rank). The driver, its termination and
jump tests, `iterations.json`, checkpoints and outputs are SDPB's. Everything is kept resident on the device: the current code
assumes that one device holds the whole iteration. Each run writes `gpu_profile.json` next to `iterations.json`:
seconds per stage and device/host high-water measurements. Options through the environment:
`SDPB_GPU_EIG_TOL` (step-length tolerance, 1e-10), `SDPB_GPU_REFINE_DY` (refinement steps in the solves with Q, 1), `SDPB_GPU_REFINE_SYSTEM` (of the corrector's Schur system: -1 automatic once the dual error is below 2^-200, 0 never, n always), `SDPB_GPU_CORRECT_XY`, `SDPB_GPU_CORRECT_SCHUR` (panel residual
correction, on), `SDPB_GPU_CHUNK_BLOCKS` (blocks per chunk, 16, which is what fits at nmax18; 0: one chunk per dimension, the fastest at nmax14, where it fits), `SDPB_GPU_STRUCTURED_P` (P by column groups, on), `SDPB_GPU_SCHUR_DIGITS` (digits of the correction's second product on Schur blocks, -1: adaptive), `SDPB_GPU_CRT_GUARD` (guard of the truncated CRT in all products, 128; 0: untruncated, bit for bit with the earlier default), `SDPB_GPU_PRODUCT_TIMING`; `SDPB_GPU_PERTURB_STEP` and
`SDPB_GPU_PERTURB_ITERATION` are the negative control of the comparison with SDPB.

### Environment settings

`--gpu`'s tuning is read from the environment of each process. The values are not echoed in `out.txt` or the
checkpoint, so record them with a run, and give every rank the same ones. `--maxMemory` and `--maxSharedMemory`
have no effect with `--gpu`: no host or device memory estimate is made.

| Variable | Default | Effect |
|---|---|---|
| `SDPB_GPU_CHUNK_BLOCKS` | 16 | blocks per chunk (0: one chunk per dimension) |
| `SDPB_GPU_SPARSE_B`, `SDPB_GPU_LEAN_P` | 0, 0 | memory modes: B kept as its nonzeros (same results); P released once Q is formed (changes the rounding) |
| `SDPB_GPU_STRUCTURED_P` | 1 | P = X B by groups of columns with the same nonzero rows |
| `SDPB_GPU_EIG_TOL` | 1e-10 | step-length tolerance |
| `SDPB_GPU_REFINE_DY`, `SDPB_GPU_REFINE_SYSTEM`, `SDPB_GPU_REFINE_X` | 1, -1, 0 | refinement: solves with Q; the corrector's Schur system (-1: once the dual error is below 2^-200); diagnostic refinement of the solves with X |
| `SDPB_GPU_CORRECT_XY`, `SDPB_GPU_CORRECT_SCHUR` | 1, 1 | panel residual correction of the factors |
| `SDPB_GPU_SCHUR_DIGITS`, `SDPB_GPU_SCHUR_LIMIT` | -1, 0 | digits of the correction's second product on Schur blocks (-1: adaptive); its limit |
| `SDPB_GPU_CRT_GUARD`, `SDPB_GPU_PQ_GUARD` | 128, 0 | guard bits of the truncated CRT (0: untruncated); the same for P = X B and Q's factor |
| `SDPB_GPU_DIGITS` | 4 L | digits per aligned operand in the products |
| `SDPB_GPU_BUDGET`, `SDPB_GPU_PANEL` | 0.7, 512 | fraction of free device memory one product may plan for; widest column panel |
| `SDPB_GPU_GEMV_THREADS`, `SDPB_GPU_GEMV_MIN_TERMS` | 262144, 16 | threads in all per matrix-vector product (0: one per output); the fewest terms a thread sums (the same result bit for bit at any split) |
| `SDPB_GPU_STRUCTURED_Q`, `SDPB_GPU_Q_BASE`, `SDPB_GPU_Q_PLANE` | 1, 2048, 2^26 | Q = P^T P by groups of rows with the same nonzero columns (0: the dense product; the same Q bit for bit); the widest diagonal block of Q formed as one full square product; on one GPU, the most entries of one product's residue plane |
| `SDPB_GPU_Q_SINGLE_PASS` | 1 | several GPUs: each rank forms P's residues once and accumulates its whole packed triangle (0: residues per tile; the same Q bit for bit, less device memory) |
| `SDPB_GPU_PRODUCT_TIMING`, `SDPB_GPU_VERBOSE`, `SDPB_GPU_CORRECTION_VERBOSE`, `SDPB_GPU_DEBUG_TIES` | unset | diagnostics: per-product timing, planner and correction messages, tie report |
| `SDPB_GPU_PERTURB_STEP`, `SDPB_GPU_PERTURB_ITERATION` | 0, 0 | test only: perturb a step, the negative control of the comparison with SDPB |
| `SDPB_GPU_CHUNK_ROWS`, `SDPB_GPU_COMM_TIMING` | 0, unset | cap on a chunk's rows (blocks times P), cutting only a dimension's largest chunks; seconds in collectives per stage in `gpu_profile.json` |
| `SDPB_GPU_CORRUPT_LEAF_ITERATION`, `SDPB_GPU_CORRUPT_LEAF_LIMBS` | 0, 0 | test only: at that iteration, factor with the low limbs of block 0 cleared (a negative control) |

The switches for Q above (`SDPB_GPU_STRUCTURED_Q`, `SDPB_GPU_Q_SINGLE_PASS`, `SDPB_GPU_GEMV_MIN_TERMS`) don't change a bit of the result. To run the code path from before them, set `SDPB_GPU_STRUCTURED_Q=0 SDPB_GPU_Q_SINGLE_PASS=0 SDPB_GPU_GEMV_MIN_TERMS=256 SDPB_GPU_Q_BASE=256`.

## Experimental distributed Q

Build with `--gpu-nccl-dir=<NCCL prefix>` (CMake: `-DSDPB_GPU_NCCL_DIR=<NCCL prefix>`) in addition to the GPU options above.
Only the NCCL headers are needed at build time. At run time set
`SDPB_GPU_NCCL_LIBRARY=/path/to/libnccl.so.2` if it is not on the library search
path. The single-rank path neither loads NCCL nor changes its factor/solve
algorithm.

Launch one MPI process per allocated GPU, for example
`mpirun -n 2 --bind-to none sdpb --gpu <the same solver options>`, and across
nodes `mpirun -n 4 --map-by ppr:2:node --bind-to none ...`, filling each node
before the next. A rank takes the GPU given by its place
among the ranks on its node; if Slurm exposes one device per rank, it uses that.
Duplicate device UUIDs are rejected. Checks cover Q's assembly, factor and
solves across two nodes; whole calls across nodes have run at nmax14 only. This
prototype runs at every L in `rns/Limb_Sizes.hxx` and requires at least
one whole chunk on each rank. Keep the same `SDPB_GPU_CHUNK_BLOCKS` and
`SDPB_GPU_CHUNK_ROWS` when comparing rank counts. `SDPB_GPU_CHUNK_ROWS` (0, off)
caps a chunk's rows, blocks times P, which cuts only the largest chunks of a
dimension and lets whole chunks spread over more ranks; it changes the chunks,
so rounding differs from the same run without it. `SDPB_GPU_COMM_TIMING=1` adds
seconds in collectives per stage to `gpu_profile.json`. Check
that a job's step really sees all its GPUs: some clusters expose only half of a node's GPUs to a step.

Known problem, not yet understood: on two nodes (2+2 ranks) `mpirun` has hung after `MPI_Abort` (three times of
three; on 1+1 it exits). Every error of a multi-node `--gpu` run (out of memory on one rank, a refused checkpoint,
fewer visible GPUs than ranks on a node, ...) ends in `MPI_Abort` from `main`, so such a job can sit idle until its
time limit. Give multi-node GPU jobs a tight time limit or a watchdog on the output.
Every rank must see the same `SDPB_GPU_*` settings (pass them with the launcher, e.g. `mpirun -x`); sdpb
compares them across ranks at start-up and refuses a difference. The block placement is written to `gpu_layout`
in the checkpoint directory, and a checkpoint is loaded only with the placement it was written with (a checkpoint
without `gpu_layout`, e.g. from a CPU run, is refused on several GPU ranks).

Whole B/P chunks and their Schur/X/Y blocks have deterministic owners.
`linalg/Distributed_Q` assigns lower Q tiles cyclically to ranks, reduces
integer residue sums before owner CRT, and retains the original owner tiles
for refinement. By default each rank accumulates the residues of its whole
packed lower triangle in one pass over its rows (`Mpf_Gram::add_structured`:
rows that share their nonzero columns form one product over those columns)
and the tiles for the collective are cut from it; this holds
primes x N (N + 1) bytes and the residue planes on every rank during the
assembly (about 15 GB at N = 7,990). `SDPB_GPU_Q_SINGLE_PASS=0` forms the
residues tile by tile instead, with less memory and about N / tile + 1 times
the residue work; the assembled Q is the same bit for bit. Corrected tiled Cholesky keeps its lower factor and diagonal
inverse tiles on their owners. It does not replicate a full Q or inverse.
Subsequent vector sums use floating arithmetic and may round differently
with the rank count. All ranks participate in the existing driver, including
checkpoints and termination; multi-rank continuation has not yet been tested.

Options: `SDPB_GPU_Q_TILE` (512), `SDPB_GPU_Q_ROW_PANEL` (8192), and the existing
`SDPB_GPU_REFINE_DY` (1). `SDPB_GPU_Q_DIGEST_PREFIX` optionally writes first-
iteration original-Q mantissa/exponent hashes to `<prefix>.rankR.txt`, with
no whole-Q host copy. Use the same digest tile size to compare rank counts.
Profiles are `gpu_profile.json` and `gpu_profile.json.rankR`. Collective
counters record logical received application payload, not measured link
traffic; reduce-scatter counts one output segment, while a ring moves G-1
segments. Host RSS is a process high water, not stage-local live memory.

Current evidence (computed once, D): nmax14/18 eight-iteration G=1/2/4 calls
have identical assembled Q tile hashes. Four-GPU device peaks are 35.0/50.7
GB, versus 82.1/146.6 GB on one GPU in the same cap16 configuration. Matched
guard128 frozen factors pass independent MPFR forward and exact GMP sampled
residual checks at every G, with corrupted factors rejected. The original
tiled full solve fails the no-worse-than-monolithic per-column gate on 2/32
and 15/32 columns. The strict relative iteration-field gate also fails for
`p-err`; other fields agree to at least 175/153 decimal digits over eight
iterations. A complete n14 G=2 call with the original residual reaches the
same dual-feasible jump at iteration145 as the original G=1 and CPU runs;
all non-p-err iteration-log fields agree to at least138/13.4 digits with those references.
Terminal out.txt scalars are separate: dualError is 1.79e-143 (G2),
2.82e-143 (G1),4.80e-168 (CPU), also failing the relative gate.
No failed gate is waived by that verdict.

The same L28 implementation also completes eight nmax22 R25+infinity
iterations on G2/G4, with matching 136 Q tile hashes, 95.3/65.6 GB device
peaks and 12.98/8.40 s/iteration. This is a resident configuration; the
one-GPU implementation with sparse B also fits that input.
The full G4 call reaches the CPU reference's dual-feasible jump at
iteration 166: 8.285 s/iteration, 65.563 GB maximum device memory, all
non-p-err iteration fields agreeing to at least 13.502 digits. Its first
eight records and all initial Q digests reproduce the prefix exactly.
The strict p-err gate still fails, and terminal dualError separately fails
(4.89e-138 GPU versus 2.31e-152 CPU). This CPU reference is a run with fewer poles
(maxNumPoles 25 and a node at infinity), not the production settings; no production CPU speedup is claimed.
One 1024-tile comparison at nmax14/18 reduces distributed iteration time
by 7% and 14–17%, with larger isolated workspaces and more collective
payload. Its forward/factor checks pass, but its full per-column gate
fails on 3/32 and 16/32 columns. The default tile remains 512. These are
single-run measurements, not production acceptance or linear scaling.

`SDPB_GPU_Q_WIDE_RESIDUAL=1` enables an optional experimental residual for
distributed Q only (default OFF; G=1 solver arithmetic stays unchanged).
It multiplies full L-limb inputs, accumulates guarded signed integers across
owners, subtracts Qv from the RHS, then rounds once. In one bounded frozen-Q test
at nmax14/18 and G=1/2, it passes exact GMP comparisons at every residual
row and the unchanged full-solve per-column gate for all32 columns at both
MPFR reference precisions. Worst full-solve errors are about 2^-894.3;
corrupted residual and factor controls fail. This is not a complete-solver
validation of the optional path, and the original residual's failure remains
recorded separately. Harness options are `--wide-residual 1 --refine 1`;
`--mode residual` exports b-Qv for an external exact checker.

Standalone harness (MPI headers/libraries must also be available to nvcc):

```
make -C src/gpu GMP_DIR=<GMP prefix> NCCL_DIR=<NCCL prefix> ARCH=90 build/distributed_q_check
mpirun -n 2 src/gpu/build/distributed_q_check --mode assembly --N 173 --tile 64
mpirun -n 2 src/gpu/build/distributed_q_check --mode factor --dir <frozen Q directory> \
  --N 3406 --tile 512 --columns 32 --crt-guard 128 --out <output prefix>
```

The assembly fixture compares every entry with the original full Gram;
`--drop-chunk 1` is its negative control. Factor exports alone are not an
accuracy gate. Compare the exported solves with the independent MPFR LDL
references of `test/ref_solve.c` at two precisions (`compare_references`;
`make mpfr-checks`), and sample the factor's residual exactly in GMP;
`--corrupt 1` must be rejected. `--mode monolithic` provides the unchanged
single-GPU inverse-product arithmetic as a separate full-solve comparator.

## Speed and memory options of the factor and P = X B

`test/factor_profile.cu` runs each on real dumps and prints digests, stages and high-water marks. On by default and
exact (the digests of L, X and P equal the plain build's on real Schur blocks and real Q):

- `factor()` works in place on views of A, L and X (`Mpf_Batch::block`): no copies per recursion level.
- `Factor_Options::lower_update`: L21 L21^T on the lower triangle, its second operand's residues from the first's cache.
- `Factor_Options::reuse_residues`: L21 X11 reads the residues of L21^T from that cache.
- `Factor_Options::structured_solve`, `Block_Factor::solve_flat`: P = X B by groups of columns of B with the same nonzero
  rows. B is 80% exact zeros in SDPB's multi-dimensional blocks and 36% of P is exactly zero (nmax14, nmax18).
- `Gram_Options::share_workspace` or `Gram_Options::workspace`: one workspace for all products on a stream, so memory is
  that of the largest product, not the sum over shapes. Pass the same `workspace` to the engines of all batch sizes.

Off by default, not exact, checked against MPFR (see [test/q_accuracy.md](test/q_accuracy.md) before relying on them):

- `Gram_Options::crt_guard_bits`: forms only the digits of each product above the rounding noise (minus the guard).
- `Factor_Options::correction_digits = -1`: the correction's second product at the digits its size needs.
- `Factor_Options::correction_limit`: a run-time estimate decides whether panels are corrected at all. Above the limit the
  factorization restarts with every panel corrected, bit for bit the corrected default, and later calls for that shape
  start corrected until `reset_correction_policy()`.

Environment: `SDPB_GPU_BUDGET` (fraction of free memory one product may plan for), `SDPB_GPU_PANEL` (widest column panel, 512).

`Block_Factor` assumes positive definite blocks (and Schur complements) at the working precision and doesn't
detect a non-positive pivot: where El::Cholesky throws, it returns meaningless factors without an error.
So `sdpb --gpu` doesn't stop where the CPU solver stops with a non-positive-definite error (a Schur block or Q that
lost definiteness, e.g. at too low a precision late in a hard run): it continues on meaningless directions.

## Tests

The checks in `test/` are built by `src/gpu/Makefile`, or by CMake with `-DSDPB_ENABLE_GPU=ON`
(targets `gpu_checks`, `mpfr_checks` and, with NCCL, `nccl_checks`, in `<build>/src/gpu/`),
and need a GPU. They are not built by waf or by CMake's default target, and not run by `test/run_all_tests.sh` or CTest.
- Self-contained, each generating its own data and comparing with independent GMP products:
  `gram_check`, `gemm_check`, `mpf_gram_check`, `mpf_gemm_check`, `chol_tile_check`, `block_factor_check`,
  `mpf_ops_check`.
  `block_factor_check` gates the residual by default; its forward-error gate is off unless `--forward-log2` is
  given.
- Need a directory of real Schur blocks, which is not distributed: `q_check`, `real_block_check`,
  `residual_check`, `gemv_compare` and `factor_profile`, and the MPFR references `ref_solve`/`compare_references` (CPU, `make mpfr-checks`) on their
  output. `test/make_spd_fixtures.py` writes synthetic SPD fixtures in the same format (see
  [test/q_accuracy.md](test/q_accuracy.md)).
- Need several ranks with NCCL (`make nccl-checks NCCL_DIR=<prefix>`): `distributed_q_check` (Q and its
  factor on G GPUs against one) and `nccl_rates` (NCCL and MPI rates at the sizes distributed Q exchanges, checking the
  data of every transfer but all-reduce, which is only timed).

## Build (standalone, beside the CPU build)

```
module load cuda/12.9.0            # any CUDA >= 12.2
make -C src/gpu GMP_DIR=<GMP prefix> ARCH=90     # 90: H100/H200, 89: L40S
src/gpu/build/gram_check --N 200 --rows 1024 --chunks 3 --digits 108   # exact check, all entries
```

Or with CMake, from the root directory: `cmake -B build -DSDPB_ENABLE_GPU=ON <the usual options>`,
then `cmake --build build --target gpu_checks`; the checks are in `build/src/gpu/`.
