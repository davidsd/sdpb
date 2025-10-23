#pragma once

#include "sdpb_util/flint.hxx"

#include <boost/noncopyable.hpp>

#include <vector>

// adopted from FLINT, fmpz_mat/mul_blas.c
// TODO: rename arguments and add description to make the code readable
struct Fmpz_Comb : boost::noncopyable
{
  fmpz_comb_t comb{};
  fmpz_comb_temp_t comb_temp{};
  const std::vector<mp_limb_t> primes;
  size_t num_primes;
  std::vector<nmod_t> mods;
  std::vector<ulong> shifts;

  Fmpz_Comb() = delete;
  // bits: Number of bits to store matrix multiplication result
  // k: number of additions.
  // If we multiply matrices A and B such that:
  //   max|A| < 2^Abits
  //   max|B| < 2^Bbits
  // then:
  //   sign = 0 if A and B are non-negative, 1 otherwise
  //   k = A.Width() = B.Height()
  //   bits = Abits + Bbits + bits(k) + sign
  Fmpz_Comb(flint_bitcnt_t bits, slong k);
  Fmpz_Comb(flint_bitcnt_t Abits, flint_bitcnt_t Bbits, int sign, slong k);
  ~Fmpz_Comb();
};
