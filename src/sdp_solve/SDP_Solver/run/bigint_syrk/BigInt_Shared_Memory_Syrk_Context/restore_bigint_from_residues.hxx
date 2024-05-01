#pragma once

#include "../fmpz/Fmpz_Comb.hxx"
#include "../fmpz/fmpz_mul_blas_util.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

inline void
restore_bigint_from_residues(const Residue_Matrices_Window<double> &window,
                             size_t i, size_t j, Fmpz_Comb &comb,
                             std::vector<mp_limb_t> &residues_buffer_temp,
                             Fmpz_BigInt &output)
{
  int sign = 1; // means that negative values are allowed
  size_t num_primes = comb.num_primes;

  residues_buffer_temp.resize(num_primes);
  for(size_t prime_index = 0; prime_index < num_primes; prime_index++)
    {
      double d = window.residues.at(prime_index)(i, j);
      ASSERT(abs(d) <= MAX_BLAS_DP_INT);
      auto &mod = comb.mods.at(prime_index);
      residues_buffer_temp.at(prime_index)
        = double_to_uint32_t_residue(d, mod);
    }

  fmpz_multi_CRT_ui(output.value, residues_buffer_temp.data(), comb.comb,
                    comb.comb_temp, sign);
}