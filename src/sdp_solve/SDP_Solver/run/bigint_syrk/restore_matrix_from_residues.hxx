#pragma once

#include <El.hpp>
#include <flint/nmod.h>
#include "Fmpz_Comb.hxx"
#include "sdpb_util/Shared_Window_Array.hxx"
#include "fmpz_mul_blas_util.hxx"
#include "fmpz_BigFloat_convert.hxx" //TODO join with prev header

inline void
restore_matrix_from_residues(const Residue_Matrices_Window<double> &window,
                             size_t i, size_t j, Fmpz_Comb &comb,
                             fmpz_t &output)
{
  int sign = 1; // means that negative values are allowed
  size_t num_primes = comb.num_primes;

  std::vector<mp_limb_t> residues(num_primes);
  for(size_t prime_index = 0; prime_index < num_primes; prime_index++)
    {
      double d = window.residues.at(prime_index)(i, j);
      assert(abs(d) <= MAX_BLAS_DP_INT);
      auto &mod = comb.mods.at(prime_index);
      residues.at(prime_index) = double_to_uint32_t_residue(d, mod);
    }

  fmpz_multi_CRT_ui(output, residues.data(), comb.comb, comb.comb_temp, sign);
}