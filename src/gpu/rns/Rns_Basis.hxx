// GPU backend for SDPB: the residue number system (RNS) used for exact products of
// big-integer matrices on INT8 tensor cores.
//
// Based on SDPviBe by Antoine Tilloy (https://github.com/AntoineTilloy/SDPviBe),
// the prototype of this backend: residues modulo ~14-15-bit primes, products on
// INT8 tensor cores, CRT reconstruction as an INT8 product with the base-256
// digits of M/m_i (SDPviBe sdpvibe/rns.py, class RNS).  Portions adapted from
// SDPviBe, Copyright (c) 2026 Antoine Tilloy, MIT License (see
// docs/GPU_backend_credits.md).
//
// Differences from SDPviBe (this file):
//  * primes below 2^14 instead of 2^15, so that a centred residue
//    r = 128 h + l and the digit sum s = h + l all fit in int8; a modular
//    product then takes three INT8 products (Karatsuba) instead of four;
//  * inputs are balanced base-256 digits (d in [-128, 127]), so residues need
//    no sign or offset corrections;
//  * the number of primes is set by the bit bound of the products, not
//    rounded to whole limbs or to groups of 16.
#pragma once

#include <cstdint>
#include <vector>

namespace sdpb_gpu
{
  struct Rns_Basis
  {
    // primes m_i, descending from 16381; all in (2^13, 2^14)
    std::vector<int32_t> primes;
    // number of primes padded to a multiple of 16 (zero residues in the pad)
    int num_padded = 0;
    // bits of M = prod m_i
    int M_bits = 0;
    // E = number of balanced base-256 digits of the CRT sums (multiple of 4)
    int crt_digits = 0;

    // Smallest basis with M > 2^bound_bits (the caller adds its margins).
    static Rns_Basis for_bound(int bound_bits);

    // Basis for exact Gram/GEMM sums of at most `rows` products of integers
    // that are given as `digits` balanced base-256 digits each
    // (|a| < 2^(8*digits - 1)): |sum| < rows * 2^(16*digits - 2), and the CRT
    // needs M > 4 |sum| (centred, with a float estimate of the quotient).
    static Rns_Basis for_products(int digits, int64_t rows);

    // Host tables (row-major as named; uploaded by Rns_Device):
    // residue table T[i][j] = 256^j mod m_i (centred), split T = 128 T_hi + T_lo
    void residue_tables(int digits, std::vector<int8_t> &T_hi,
                        std::vector<int8_t> &T_lo) const;
    // inverses inv[i] = (M/m_i)^-1 mod m_i
    std::vector<int32_t> crt_inverses() const;
    // D[j][i] = j-th balanced base-256 digit of M/m_i, for j < crt_digits,
    // i < num_padded (zero for padding primes); and the digits of M
    void crt_tables(std::vector<int8_t> &D, std::vector<int32_t> &M_digits) const;
    // Single-product CRT: with balanced base-128
    // digits d_ij of M/m_i and y_i = 128 h_i + l_i, sum_i y_i M/m_i =
    // sum_j Z_j 128^j where Z_j = sum_i l_i d_ij + sum_i h_i d_i,j-1: one INT8
    // product of [l | h] (2 num_padded per entry) with D2, D2[j][i] = d_ij,
    // D2[j][num_padded + i] = d_i,j-1 (j < crt_digits_128()); and M's base-128
    // digits.
    int crt_digits_128() const;
    void crt_tables_128(std::vector<int8_t> &D2,
                        std::vector<int32_t> &M_digits) const;
  };
}
