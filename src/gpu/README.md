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
- chunked streaming of per-block arrays between host and device, and NCCL multi-GPU (planned).

Device-float products first round operands to the chosen per-column
alignment. RNS accumulation is exact for those aligned integers; the result
is then rounded to a floating mantissa. The real-Q accuracy/cost tradeoff and
explicit options are documented in [the regression recipe](test/q_accuracy.md).

## What is here (products and block factorization)

- `rns/Rns_Basis.{hxx,cxx}`: primes in (2^13, 2^14), residue and CRT tables (host, GMP).
- `rns/Gram_Accumulator.{hxx,cu}`: exact Q = Σ_chunks AᵀA for big-integer A (balanced base-256 digits),
  residues and products through cuBLASLt INT8 GEMMs (wgmma on Hopper), three INT8 products per prime
  (Karatsuba; SDPviBe uses four), the lower triangle cut into recursive rectangles, uint16 residues of the
  packed lower triangle kept across chunks, CRT once in `finish()`.
- `rns/Mpf_Gram.{hxx,cu}`: the same on device floats (SDPviBe's format, column-major): column maxima,
  residues straight from the 32-bit limbs, exact Gram, each entry of Q rounded once to a float.
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

`mpirun -n 1 sdpb --gpu <options>` runs the iterations on one GPU (one MPI rank). The driver, its termination and
jump tests, `iterations.json`, checkpoints and outputs are SDPB's. Everything is kept resident on the device: the current code
assumes that one device holds the whole iteration. Each run writes `gpu_profile.json` next to `iterations.json`:
seconds per stage and the device memory at the end of each stage. Options through the environment:
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
| `SDPB_GPU_GEMV_THREADS` | 262144 | threads in all per matrix-vector product (0: one per output) |
| `SDPB_GPU_PRODUCT_TIMING`, `SDPB_GPU_VERBOSE`, `SDPB_GPU_CORRECTION_VERBOSE`, `SDPB_GPU_DEBUG_TIES` | unset | diagnostics: per-product timing, planner and correction messages, tie report |
| `SDPB_GPU_PERTURB_STEP`, `SDPB_GPU_PERTURB_ITERATION` | 0, 0 | test only: perturb a step, the negative control of the comparison with SDPB |

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
(targets `gpu_checks`, `mpfr_checks`, in `<build>/src/gpu/`),
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

## Build (standalone, beside the CPU build)

```
module load cuda/12.9.0            # any CUDA >= 12.2
make -C src/gpu GMP_DIR=<GMP prefix> ARCH=90     # 90: H100/H200, 89: L40S
src/gpu/build/gram_check --N 200 --rows 1024 --chunks 3 --digits 108   # exact check, all entries
```

Or with CMake, from the root directory: `cmake -B build -DSDPB_ENABLE_GPU=ON <the usual options>`,
then `cmake --build build --target gpu_checks`; the checks are in `build/src/gpu/`.
