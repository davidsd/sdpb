#include "residues.hxx"
#include "sdpb_util/assert.hxx"

// Helper functions adopted from FLINT library, mat_mul_blas.c

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

// out[l * stride] := input mod comb.primes[l], in the range (-p/2; p/2]
void fmpz_multi_mod_uint32_stride(double *out, slong stride,
                                  const fmpz_t &input, Fmpz_Comb &comb)
{
  auto &residues = comb.residues_temp;
  fmpz_multi_mod_ui(residues.data(), input, comb.comb, comb.comb_temp);
  for(size_t l = 0; l < comb.num_primes; ++l)
    out[l * stride] = uint32_t_residue_to_double(residues[l], comb.primes[l]);
}

// Convert BigInt to and from residues

void bigint_to_residues(const Fmpz_BigInt &input, Fmpz_Comb &comb,
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