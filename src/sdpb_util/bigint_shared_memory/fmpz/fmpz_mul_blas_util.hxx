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

inline double _reduce_uint32(mp_limb_t a, nmod_t mod)
{
  mp_limb_t r;
  NMOD_RED(r, a, mod);
  return (uint32_t)r;
}

inline double _reduce_double(mp_limb_t a, nmod_t mod)
{
  return uint32_t_residue_to_double(_reduce_uint32(a, mod), mod.n);
}

inline void
fmpz_multi_mod_uint32_stride(double *out, slong stride, const fmpz_t &input,
                             const Fmpz_Comb &comb)
{
  const fmpz_comb_t &C = comb.comb;
  const fmpz_comb_temp_t &CT = comb.comb_temp;

  slong i, k, l;
  fmpz *A = CT->A;
  mod_lut_entry *lu;
  slong *offsets;
  slong klen = C->mod_klen;
  fmpz_t ttt;

  /* high level split */
  if(klen == 1)
    {
      *ttt = A[0];
      A[0] = *input;
    }
  else
    {
      _fmpz_multi_mod_precomp(A, C->mod_P, input, -1, CT->T);
    }

  offsets = C->mod_offsets;
  lu = C->mod_lu;

  for(k = 0, i = 0, l = 0; k < klen; k++)
    {
      slong j = offsets[k];

      for(; i < j; i++)
        {
          /* mid level split: depends on FMPZ_MOD_UI_CUTOFF */
          mp_limb_t t = fmpz_get_nmod(A + k, lu[i].mod);

          /* low level split: 1, 2, or 3 small primes */
          if(lu[i].mod2.n != 0)
            {
              FLINT_ASSERT(l + 3 <= C->num_primes);
              out[l * stride] = _reduce_double(t, lu[i].mod0);
              l++;
              out[l * stride] = _reduce_double(t, lu[i].mod1);
              l++;
              out[l * stride] = _reduce_double(t, lu[i].mod2);
              l++;
            }
          else if(lu[i].mod1.n != 0)
            {
              ASSERT(l + 2 <= C->num_primes);
              out[l * stride] = _reduce_double(t, lu[i].mod0);
              l++;
              out[l * stride] = _reduce_double(t, lu[i].mod1);
              l++;
            }
          else
            {
              ASSERT(l + 1 <= C->num_primes);
              out[l * stride]
                = uint32_t_residue_to_double((uint32_t)(t), comb.primes.at(l));
              l++;
            }
        }
    }

  ASSERT(l == C->num_primes);

  if(klen == 1)
    A[0] = *ttt;
}
