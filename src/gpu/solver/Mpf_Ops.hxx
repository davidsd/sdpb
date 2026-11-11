// GPU backend for SDPB: element kernels on device floats, a direct
// multiprecision matrix-vector product, reductions, and the solver-specific
// assemblies (Schur complement, constraint-matrix weighted sums, pairing
// diagonals).
//
// Based on SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe),
// the prototype of this backend: the float format and elementwise arithmetic
// follow sdpvibe/mpf.py (mpf_mul, mpf_add, total, max_abs, to_float), the
// assemblies follow sdpvibe/solver.py (cmws, schur, seg_diag).  Portions adapted
// from SDPviBe, Copyright (c) 2026 Antoine Tilloy, MIT License
// (docs/GPU_backend_credits.md).
//
// New here: all the assemblies are single launches over the batch (one launch
// per segment pair in SDPviBe), and the matrix-vector products are schoolbook
// products on CUDA cores, without residues: a product with one vector needs
// the residues of the whole matrix in residue arithmetic, which costs more than
// the schoolbook products themselves.
#pragma once

#include "Host_Float.hxx"
#include "../linalg/Block_Factor.hxx"

#include <cstdint>
#include <vector>

namespace sdpb_gpu
{
  // O[i] = A[i] B[i]; with B 1 x 1 x 1, O[i] = A[i] B[0]
  void dev_mul(const Mpf_Batch &A, const Mpf_Batch &B, Mpf_Batch &O,
               cudaStream_t s);
  // O = -A (O may be A)
  void dev_neg(const Mpf_Batch &A, Mpf_Batch &O, cudaStream_t s);
  // A *= 2^k, exactly
  void mpf_shift_exp(Mpf_Batch &A, int64_t k, cudaStream_t s);
  // M(i, j) = M(j, i) for i < j
  void mpf_lower_to_upper(Mpf_Batch &M, cudaStream_t s);
  // zero the entries of the n x n blocks (n = dim h) outside the valid rows and
  // columns: with block size h, only (i % h) < hj[z] and (j % h) < hj[z] are valid
  void mpf_mask_blocks(Mpf_Batch &M, int h, const int *d_hj, cudaStream_t s);
  // M(i, i) += a on the valid rows (a 1 x 1 x 1)
  void mpf_add_diag(Mpf_Batch &M, const Mpf_Batch &a, int h, const int *d_hj,
                    cudaStream_t s);
  // dst (npad x npad) = [src 0; 0 I], src n x n
  void mpf_pad_identity(const Mpf_Batch &src, Mpf_Batch &dst, cudaStream_t s);
  // D (n x 1) = diagonal of M (n x n)
  void mpf_diagonal(const Mpf_Batch &M, Mpf_Batch &D, cudaStream_t s);
  // Vx (batch, dim h x dim K): V (batch, h x K) on the diagonal, blocks of size h x K
  void mpf_kron_identity(const Mpf_Batch &V, Mpf_Batch &Vx, int dim,
                         cudaStream_t s);
  // R (batch * copies, rows x cols) = copies of A (batch, rows x cols):
  // R(b + batch * c) = A(b)
  void mpf_repeat(const Mpf_Batch &A, Mpf_Batch &R, cudaStream_t s);

  // B stored as its nonzero entries, after Zechuan Zheng's SparseMPF format:
  // pos[i] is the flat index of nonzero i, ascending. An entry is packed
  // unless every bit of it says zero (exponent MPF_EZ and mantissa 0), so that
  // unpacking gives back the same bytes whatever else an entry holds.
  void mpf_pack_sparse(const int32_t *m, const int64_t *e, int64_t entries, int L,
                       std::vector<uint32_t> &pos, std::vector<int32_t> &vm,
                       std::vector<int64_t> &ve);
  // A = exact zeros, except A(pos[i]) = v(i); v is (1, nnz x 1)
  void mpf_unpack_sparse(const Mpf_Batch &v, const uint32_t *pos, Mpf_Batch &A,
                         cudaStream_t s);

  // Segments of a block of matrix dimension dim: s = c (c + 1) / 2 + r, r <= c.
  // v (nb, P x 1), P = nseg K: v(q, s K + k) = sum over parity of
  // M(parity nb + q)(r K + k, c K + k).  M: (2 nb, dim K x dim K)
  void mpf_seg_diag(const Mpf_Batch &M, Mpf_Batch &v, int K, int dim,
                    cudaStream_t s);
  // VD (2 nb nseg, K x h): VD(z + 2 nb s)(k, i) = a(q, s K + k) VT(z)(k, i),
  // z = parity nb + q.  VT: (2 nb, K x h); a: (nb, P x 1)
  void mpf_scale_segments(const Mpf_Batch &VT, const Mpf_Batch &a, Mpf_Batch &VD,
                          int nseg, cudaStream_t s);
  // R (2 nb, dim h x dim h) from C (2 nb nseg, h x h): the block (r, c), r <= c,
  // is C(z + 2 nb s), times 1/2 if r != c, and the block (c, r) its transpose
  void mpf_cmws_assemble(const Mpf_Batch &C, Mpf_Batch &R, int h, int dim,
                         cudaStream_t s);
  // S (nb, P x P), P = nseg K (SDPB's compute_schur_complement): AXI, AY are
  // (2 nb, dim K x dim K).  Padded points (k >= Kj[q]) get a unit diagonal.
  void mpf_schur(const Mpf_Batch &AXI, const Mpf_Batch &AY, Mpf_Batch &S, int K,
                 int dim, const int *d_Kj, cudaStream_t s);

  // out[b ob + o oo] = sum_{tb < ntb, ti < nti} A[ia] v[iv], for b < batch and
  // o < no, where ia = b ab + o ao + tb atb + ti ati and
  // iv = b vb + o vo + tb vtb + ti vti.  A and v are device floats in limb
  // planes; out's plane is oplane.  Rounded once per entry (error about 2^-NB
  // relative to the largest term of each sum).
  struct Gemv_Spec
  {
    int batch = 1;
    int64_t no = 0;
    int64_t ab = 0, ao = 0, atb = 0, ati = 0;
    int64_t vb = 0, vo = 0, vtb = 0, vti = 0;
    int64_t ntb = 1, nti = 0;
    int64_t ob = 0, oo = 1;
  };
  // Threads per output: about gemv_thread_target / outputs, at most 64; 0: one
  // thread per output.
  extern int gemv_thread_target;
  void mpf_gemv(const Mpf_Batch &A, const Mpf_Batch &v, Mpf_Batch &out,
                const Gemv_Spec &spec, cudaStream_t s);

  // Reductions (they synchronize the stream).  max |A[i]| and the index of an
  // entry attaining it (zero if all entries are zero).
  Host_Float mpf_max_abs(const Mpf_Batch &A, cudaStream_t s);
  // sum of all the entries (exact to about 2^-(NB + 32) of the largest entry);
  // the mantissa may have more limbs than a float
  Host_Float mpf_total(const Mpf_Batch &A, cudaStream_t s);

  // Eigenvalue input: for each block b, the entries as doubles scaled by
  // 2^-scale[b], with scale[b] the block's largest e + NB (0 if the block is
  // zero).  D: batch * rows * cols doubles, in the batch layout.
  void mpf_to_double(const Mpf_Batch &A, double *d_out, int64_t *scale,
                     cudaStream_t s);
  // limb planes and exponents of all the entries, device to host
  void mpf_download(const Mpf_Batch &A, int32_t *m, int64_t *e, cudaStream_t s);
  void mpf_upload(const int32_t *m, const int64_t *e, Mpf_Batch &A,
                  cudaStream_t s);
}
