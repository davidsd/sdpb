#pragma once

#include "Fmpz_BigInt.hxx"
#include "Fmpz_Comb.hxx"

#include <vector>

void bigint_to_residues(const Fmpz_BigInt &input, const Fmpz_Comb &comb,
                        double *residues, size_t stride);

void residues_to_bigint(const double *residues, size_t stride, Fmpz_Comb &comb,
                        std::vector<mp_limb_t> &residues_buffer_temp,
                        Fmpz_BigInt &output);