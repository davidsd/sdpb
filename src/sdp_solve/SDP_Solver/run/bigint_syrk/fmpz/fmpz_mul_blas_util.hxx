#pragma once

#include "Fmpz_Comb.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

// Code adopted from FLINT library, mat_mul_blas.c

#define MAX_BLAS_DP_INT (UWORD(1) << 53)

// see mul_blas.c, _tod_worker()
// n is the prime
// input value is in the range [0;n)
// output is in the range (-n/2; n/2]
inline double uint32_t_residue_to_double(uint32_t value, mp_limb_t n)
{
  // return (int32_t)(value - (n & (-(uint32_t)((int32_t)(n/2 - value) < 0))));

  // This is equivalent to the above formula, but more readable:
  if(n / 2 < value)
    return (int32_t)(value - (n & UINT32_MAX));
  else
    return value;
}

inline uint32_t double_to_uint32_t_residue(double value, nmod_t mod)
{
  ulong shift = ((2 * MAX_BLAS_DP_INT) / mod.n) * mod.n;
  mp_limb_t r;
  slong a = (slong)value;
  mp_limb_t b = (a < 0) ? a + shift : a;
  NMOD_RED(r, b, mod); // r := b % mod.n
  return (uint32_t)r;
}

// out[l * stride] := input mod comb.primes[l], in the range (-p/2; p/2]
inline void fmpz_multi_mod_uint32_stride(double *out, slong stride,
                                         const fmpz_t &input, Fmpz_Comb &comb)
{
  auto &residues = comb.residues_temp;
  fmpz_multi_mod_ui(residues.data(), input, comb.comb, comb.comb_temp);
  for(size_t l = 0; l < comb.num_primes; ++l)
    out[l * stride] = uint32_t_residue_to_double(residues[l], comb.primes[l]);
}
