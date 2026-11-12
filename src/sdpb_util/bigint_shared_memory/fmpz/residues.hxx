#pragma once

#include "Fmpz_BigInt.hxx"
#include "Fmpz_Comb.hxx"

#include <vector>

inline constexpr ulong MAX_BLAS_DP_INT = UWORD(1) << 53;

// A residue r in [0, n) as a double in (-n/2, n/2] (see FLINT's mul_blas.c)
double uint32_t_residue_to_double(uint32_t value, mp_limb_t n);
// An integer-valued double |value| <= MAX_BLAS_DP_INT reduced mod mod.n,
// in [0, mod.n)
uint32_t double_to_uint32_t_residue(double value, nmod_t mod);

void bigint_to_residues(const Fmpz_BigInt &input, Fmpz_Comb &comb,
                        double *residues, size_t stride);

void residues_to_bigint(const double *residues, size_t stride, Fmpz_Comb &comb,
                        std::vector<mp_limb_t> &residues_buffer_temp,
                        Fmpz_BigInt &output);