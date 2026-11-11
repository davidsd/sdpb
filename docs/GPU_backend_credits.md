# GPU backend: what comes from SDPviBe

The GPU backend is based on SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe, MIT
License), the prototype it follows. This file lists what was taken from it and where, so that releases and
papers can credit it accurately. Code that adapts SDPviBe's code (not only its ideas) keeps its notice:

    Portions adapted from SDPviBe, Copyright (c) 2026 Antoine Tilloy, MIT License.

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
    associated documentation files (the "Software"), to deal in the Software without restriction,
    including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions: The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
    WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    OTHER DEALINGS IN THE SOFTWARE.

SDPviBe was read at commit 0d02884735158b78e3a9ddd16c070ab20de90e63.

| Backend file | From SDPviBe (file: what) | Changed here |
|---|---|---|
| `src/gpu/rns/Rns_Basis.{hxx,cxx}` | `sdpvibe/rns.py`, class `RNS`: residues modulo small primes, tables of 256^j mod m_i split into int8 halves, CRT tables of the base-256 digits of M/m_i, quotient k from Σ y_i/m_i | primes below 2^14 (so h + l fits in int8: three products per prime), balanced-digit inputs (no sign/offset constants), prime count from the bit bound |
| `src/gpu/rns/Gram_Accumulator.{hxx,cu}` | `sdpvibe/rns.py` (`_residue_kernel`, `_crt_dot_kernel`, `crt_carry`), `sdpvibe/batched.py` (`_syrk_lower_kernel`), `sdpvibe/mpf.py` (`GramAcc`: residues of every chunk's PᵀP accumulated, one CRT) | Karatsuba (3 INT8 products instead of 4), cuBLASLt INT8 GEMMs instead of Triton kernels, recursive lower-triangle rectangles, packed uint16 residue accumulator |
| `src/gpu/rns/rns_kernels.cuh` | `sdpvibe/rns.py` (`_residue_kernel`, `_crt_dot_kernel`, `crt_carry`, `crt_pre`) | CUDA C++ kernels around cuBLASLt INT8 GEMMs; balanced digits; Karatsuba split |
| `src/gpu/rns/Mpf_Gram.{hxx,cu}` | `sdpvibe/mpf.py` (the MPF float format: L 32-bit limbs and an int64 exponent, `normalize`; `GramAcc`'s column alignment), `sdpvibe/limbops.py` (`shift_round`) | exact column maxima from a first pass over resident chunks instead of running references with 28 bits of headroom and flushes |
| `src/gpu/rns/Rns_Gemm.{hxx,cu}` | `sdpvibe/batched.py` (`bgemm_nt`, `_bgemm_kernel`: products batched over equal-size blocks) | Karatsuba, cuBLASLt, column panels with one CRT per panel |
| `src/gpu/linalg/basecase.cuh` | `sdpvibe/basecase.py` (`_SRC`: `mulfx`, `rsqrt_fx`, warp-cooperative `wmul`/`wrsqrt`, `chol_tiles`, `trinv_tiles`): the CUDA C code itself | NVRTC macros (WL, LF, BT) turned into C++ templates; arithmetic unchanged |
| `src/gpu/linalg/Block_Factor.{hxx,cu}` | `sdpvibe/linalg.py` (`chol_inv`: L21 = A21 X11ᵀ, S = A22 − L21 L21ᵀ, X21 = −X22 L21 X11), `sdpvibe/mpf.py` (`Chol`: power-of-two equilibration; `mpf_add`, `shift_mag`), `sdpvibe/basecase.py` leaves | exact batched integer products after float alignment (Mpf_Gemm); P = X B with refinement; whole-matrix scaling and panel residual correction; optional adaptive correction precision. The tested `trsm_right` ordering is not included. |
| `src/gpu/rns/rns_kernels.cuh` (`residues_fused`, `crt_fused`) | the fused residue and CRT kernels of `sdpvibe/rns.py` (`_residue_kernel`, `_crt_dot_kernel` + `crt_carry`), as structure | CUDA WMMA instead of Triton; finish and carries in shared memory and registers |
| `src/gpu/linalg/leaf_fused.cuh` | `sdpvibe/basecase.py`: the fixed-point format (WL limbs, LF fractional limbs), the truncation of products below column LF − 2, the precision-doubling Newton schedule for the pivot's reciprocal square root from a double, and `Tile_Kernels` (`mulfx`, `to_double`, `from_double`) used as ported | the Cholesky and the inverse fused into one right-looking pass (forward substitution for LX = I alongside the factorization), products spread over lanes with shuffle-resolved carries, a fixed fully unrolled product for operands ≤ 1, look-ahead pivot |
| `src/gpu/rns/Rns_Basis.cxx` (`crt_tables_128`), `rns_kernels.cuh` (`crt_carry_128`), `Rns_Gemm.cu` (`gemm_reduce_y`) | the CRT of `sdpvibe/rns.py` (`crt_pre`, `_crt_dot_kernel`, `crt_carry`): y_i, the quotient k, x = Σ y_i M/m_i − kM | base-128 digit tables so that the CRT is one INT8 product; the CRT inputs formed in the modular reduction; the quotient accumulated in single precision |
| `src/gpu/bench/sdpvibe_gram_bench.py`, `sdpvibe_block_bench.py`, `sdpvibe_real_bench.py`, `sdpvibe_q_bench.py` | run SDPviBe's own `GramAcc`, `gemm`, `Chol`, `solve_panels` unchanged, for matched timing (real block/Q dumps where named) | — |

| `src/gpu/solver/Mpf_Ops.{hxx,cu}` | `sdpvibe/mpf.py` (`mpf_mul`, `mpf_add`, `normalize`, `total`, `max_abs`, `to_float`), `sdpvibe/solver.py` (`cmws`, `schur`, `seg_diag`, `masked`, `add_diag`) | templates with the limbs in registers; assemblies as one launch over the batch; matrix-vector products as schoolbook sums on CUDA cores (new) |
| `src/gpu/solver/Gpu_Solver.{hxx,cu}` | `sdpvibe/solver.py` (`Solver.run`: the passes and their order; `SDP._build`: blocks batched by dimension with inert padding; `min_eig`: float64 step lengths, the error bound (n + 2) eps ||M||_F, the tolerance) | explicit inverses instead of triangular solves; a block the float64 eigenvalue cannot resolve goes to SDPB's own multiprecision eigenvalue on the host, not to SDPviBe's `refine_min_eig` |

Planned (to be listed here as they are ported): the tiled Cholesky of Q (`blockchol.py`), the fixed-point
recursive Cholesky with triangular solves and leaf refinement (`linalg.chol`, `batched.py`), `refine_min_eig`, the chunk/streaming planner
(`solver.py` `Layout`, `store.py`), multi-GPU over NCCL (`dist.py`).
