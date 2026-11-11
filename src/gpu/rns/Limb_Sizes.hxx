// GPU backend for SDPB: the float sizes that have kernels compiled for them.
//
// A float has L limbs of 32 bits.  The fastest kernels keep every limb in a
// register, so L is a template argument of theirs.  This list is the only place
// that names the sizes: even L from 24 to 40. By the rule in Gpu_Compute
// (L = ceil((p + 2) / 32) + 1, rounded up to a listed size) that accepts every
// --precision p up to 1246 bits; below 704 bits L is rounded up to 24.  To cover
// more, add to the list: nothing else depends on the values.  Every entry costs
// compile time, in Mpf_Ops.cu and Block_Factor.cu above all.
//
// An L not in the list still works in Block_Factor and Mpf_Gemm: the leaf runs
// at the next compiled size, with more fractional bits than asked for, and
// residues and rounding use the kernels whose L is a run-time argument (slower,
// the same bits).  Mpf_Ops has no such path, so the solver rounds L up to a
// compiled size: rounded_L below.
#pragma once

// F(L, ...) for every compiled L
#define SDPB_GPU_EACH_L(F, ...)                                               \
  F(24, __VA_ARGS__)                                                          \
  F(26, __VA_ARGS__)                                                          \
  F(28, __VA_ARGS__)                                                          \
  F(30, __VA_ARGS__)                                                          \
  F(32, __VA_ARGS__)                                                          \
  F(34, __VA_ARGS__)                                                          \
  F(36, __VA_ARGS__)                                                          \
  F(38, __VA_ARGS__)                                                          \
  F(40, __VA_ARGS__)

namespace sdpb_gpu
{
  namespace limb_sizes_detail
  {
#define SDPB_GPU_LIMB_ENTRY_(L, ...) L,
    constexpr int list[] = {SDPB_GPU_EACH_L(SDPB_GPU_LIMB_ENTRY_)};
#undef SDPB_GPU_LIMB_ENTRY_
  }

  constexpr int largest_compiled_L()
  {
    int best = 0;
    for(int v : limb_sizes_detail::list)
      best = v > best ? v : best;
    return best;
  }

  // the smallest compiled size >= L; 0 if there is none
  constexpr int rounded_L(int L)
  {
    int best = 0;
    for(int v : limb_sizes_detail::list)
      if(v >= L && (best == 0 || v < best))
        best = v;
    return best;
  }

  constexpr bool is_compiled_L(int L) { return rounded_L(L) == L; }

  // The leaves' wide fixed point: LF fractional limbs and at least 24 integer
  // limbs, a multiple of four in all (the leaf moves limbs four at a time)
  constexpr int leaf_limbs(int LF) { return (LF + 24 + 3) / 4 * 4; }

  // Register widths of the rounding kernel for a product of W limbs, rounded to
  // L limbs: it must hold L + 1 limbs, and the largest a product of 4 L digits
  // can have (2 L plus the CRT's headroom; see Rns_Basis::for_bound).
  constexpr int rounding_regs_narrow(int L)
  {
    const int r = (L + 8) / 8 * 8;
    return r > 32 ? r : 32;
  }
  constexpr int rounding_regs_wide(int L)
  {
    const int r = (2 * L + 11) / 8 * 8;
    return r > 64 ? r : 64;
  }
}
