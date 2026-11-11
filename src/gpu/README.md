# GPU backend for SDPB (in development)

**The GPU backend is based on [SDPviBe](https://github.com/AntoineTilloy/SDPviBe) by Antoine Tilloy, the
prototype it follows.** SDPviBe re-implements SDPB 3.1's interior-point method on GPUs; this backend brings its
algorithms into SDPB itself, in C++/CUDA, so that SDPB stays one codebase. Taken from SDPviBe (details and file
mapping in [docs/GPU_backend_credits.md](../../docs/GPU_backend_credits.md)):

- exact multiprecision matrix products in a residue number system (RNS), with the products on INT8 tensor
  cores and the reconstruction (CRT) as an INT8 product with the digits of M/m_i;
- the Gram matrix Q = PᵀP accumulated in residue form over chunks, with one reconstruction per iteration;
- batched, equilibrated fixed-point Cholesky and triangular inverses over SDP blocks of equal size;
- step lengths from float64 eigenvalues with an error bound and multiprecision refinement (planned);
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

`Block_Factor` assumes positive definite blocks (and Schur complements) at the working precision and doesn't
detect a non-positive pivot: where El::Cholesky throws, it returns meaningless factors without an error.

## Tests

The checks in `test/` are built by `src/gpu/Makefile`, or by CMake with `-DSDPB_ENABLE_GPU=ON`
(targets `gpu_checks`, `mpfr_checks`, in `<build>/src/gpu/`),
and need a GPU. They are not built by waf or by CMake's default target, and not run by `test/run_all_tests.sh` or CTest.
- Self-contained, each generating its own data and comparing with independent GMP products:
  `gram_check`, `gemm_check`, `mpf_gram_check`, `mpf_gemm_check`, `chol_tile_check`, `block_factor_check`.
  `block_factor_check` gates the residual by default; its forward-error gate is off unless `--forward-log2` is
  given.
- Need a directory of real Schur blocks, which is not distributed: `q_check`, `real_block_check` and
  `residual_check`, and the MPFR references `ref_solve`/`compare_references` (CPU, `make mpfr-checks`) on their
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
