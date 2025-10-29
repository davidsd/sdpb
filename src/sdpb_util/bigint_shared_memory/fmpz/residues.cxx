#include "residues.hxx"
#include "sdpb_util/assert.hxx"

// Helper functions adopted from FLINT library, mat_mul_blas.c

#define MAX_BLAS_DP_INT (UWORD(1) << 53)

// see mul_blas.c, _tod_worker()
// n is the prime
// input value is in the range [0;n)
// output is in the range (-n/2; n/2]
double uint32_t_residue_to_double(const uint32_t value, const mp_limb_t n)
{
  // return (int32_t)(value - (n & (-(uint32_t)((int32_t)(n/2 - value) < 0))));

  // This is equivalent to the above formula, but more readable:
  if(n / 2 < value)
    return (int32_t)(value - (n & UINT32_MAX));
  else
    return value;
}

uint32_t double_to_uint32_t_residue(double value, nmod_t mod)
{
  ulong shift = ((2 * MAX_BLAS_DP_INT) / mod.n) * mod.n;
  mp_limb_t r;
  slong a = (slong)value;
  mp_limb_t b = (a < 0) ? a + shift : a;
  NMOD_RED(r, b, mod); // r := b % mod.n
  return (uint32_t)r;
}

double _reduce_uint32(mp_limb_t a, nmod_t mod)
{
  mp_limb_t r;
  NMOD_RED(r, a, mod);
  return (uint32_t)r;
}

double _reduce_double(const mp_limb_t a, nmod_t mod)
{
  return uint32_t_residue_to_double(_reduce_uint32(a, mod), mod.n);
}

void fmpz_multi_mod_uint32_stride(double *out, slong stride,
                                  const fmpz_t &input, const Fmpz_Comb &comb)
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

// Convert BigInt to and from residues

void bigint_to_residues(const Fmpz_BigInt &input, const Fmpz_Comb &comb,
                        double *residues, const size_t stride)
{
  fmpz_multi_mod_uint32_stride(residues, stride, input.value, comb);
}

void residues_to_bigint(const double *residues, const size_t stride,
                        Fmpz_Comb &comb,
                        std::vector<mp_limb_t> &residues_buffer_temp,
                        Fmpz_BigInt &output)
{
  constexpr int sign = 1; // means that negative values are allowed
  const size_t num_primes = comb.num_primes;

  residues_buffer_temp.resize(num_primes);
  for(size_t prime_index = 0; prime_index < num_primes; prime_index++)
    {
      const double d = residues[prime_index * stride];
      ASSERT(abs(d) <= MAX_BLAS_DP_INT);
      const auto &mod = comb.mods.at(prime_index);
      residues_buffer_temp.at(prime_index)
        = double_to_uint32_t_residue(d, mod);
    }
  if(fmpz_cmpabs(output.value, comb.primes_product.value) >= 0)
    {
      El::BigFloat output_value;
      El::BigFloat primes_product_value;
      output.to_BigFloat(output_value);
      comb.primes_product.to_BigFloat(primes_product_value);
      RUNTIME_ERROR(
        "restore_from_residues overflow: restored value = ", output_value,
        " exceeds product of comb primes = ", primes_product_value);
    }
  fmpz_multi_CRT_ui(output.value, residues_buffer_temp.data(), comb.comb,
                    comb.comb_temp, sign);
}