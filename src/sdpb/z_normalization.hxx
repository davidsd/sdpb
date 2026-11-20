#pragma once

// z.txt: the PMP functional z = y with z_k inserted at
// k = max_normalization_index(n), where z_k is chosen so that n.z == 1.
//
// n.z == 1 must hold for the written decimals, not only in memory.
// z_k = (1 - sum_{i!=k} n_i z_i) / n_k cancels terms of size max|n_i z_i|
// down to O(1), so forming it at SDPB's precision P loses
// log2 max|n_i z_i| bits. After a dual feasible jump |z| can exceed 2^P:
// at P = 128 a z.txt with |z| ~ 1e84 had n.z = -9.4e52.
//
// So z_k is formed from the decimals written for the other components (the
// strings of y.txt) and from n (SDP::normalization, read from
// normalization.json at the precision of its decimals), at
//   P' = P + log2(1 + max|n_i z_i|) + log2(N) + 64 bits,
// and written with all P' bits. Each rounding is then <= 2^-P' relative to
// a term of size <= max|n_i z_i|, so the written decimals satisfy
// |n.z - 1| <~ 2^-(P + 60).

#include "pmp/max_normalization_index.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Boost_Float.hxx"

#include <El.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace z_normalization
{
  // y's entries as El::Print writes them to y.txt
  inline std::vector<std::string>
  printed_entries(const El::Matrix<El::BigFloat> &y)
  {
    std::ostringstream os;
    El::Print(y, "", "\n", os);
    std::istringstream is(os.str());
    std::vector<std::string> entries{std::istream_iterator<std::string>(is),
                                     std::istream_iterator<std::string>()};
    ASSERT_EQUAL(entries.size(), static_cast<size_t>(y.Height()));
    return entries;
  }

  // Precision in bits for computing n.z of the decimals z (z_k excluded):
  // P' = P + log2(1 + max_{i!=k} |n_i z_i|) + log2(N) + 64
  inline mp_bitcnt_t
  working_precision(const std::vector<std::string> &z,
                    const std::vector<Boost_Float> &normalization,
                    const size_t k, const mp_bitcnt_t precision)
  {
    int max_exponent = 0;
    for(size_t i = 0; i < z.size(); ++i)
      if(i != k)
        {
          int exponent = 0;
          frexp(Boost_Float(normalization[i] * Boost_Float(z[i])), &exponent);
          max_exponent = std::max(max_exponent, exponent);
        }
    return precision + max_exponent + 1
           + static_cast<mp_bitcnt_t>(std::ceil(std::log2(z.size()))) + 64;
  }

  // The decimals of z, given the decimals of y, the normalization n,
  // the eliminated index k, and SDPB's precision in bits.
  inline std::vector<std::string>
  z_entries(const std::vector<std::string> &y,
            const std::vector<Boost_Float> &normalization, const size_t k,
            const mp_bitcnt_t precision)
  {
    const size_t N = normalization.size();
    ASSERT_EQUAL(N, y.size() + 1);
    ASSERT(k < N, DEBUG_STRING(k), DEBUG_STRING(N));
    ASSERT(normalization[k] != 0, "normalization[", k, "] == 0");

    std::vector<std::string> z(y);
    z.insert(z.begin() + k, "");

    Scoped_Boost_Float_Precision scoped_precision(
      working_precision(z, normalization, k, precision));
    // z_k = (1 - sum_{i!=k} n_i z_i) / n_k
    Boost_Float z_k = 1;
    for(size_t i = 0; i < N; ++i)
      if(i != k)
        z_k -= normalization[i] * Boost_Float(z[i]);
    z_k /= normalization[k];
    // all digits of z_k, in the format of the other entries (as El::Print
    // writes a BigFloat: gmpxx's operator<< with default flags)
    z[k] = z_k.str(0, std::ios_base::fmtflags());
    return z;
  }

  // |n.z - 1| for the decimals z (z_k included), at the working precision.
  // A self-check of z_entries(): it catches decimals that do not read back
  // as written, not an error in n.
  inline Boost_Float
  normalization_error(const std::vector<std::string> &z,
                      const std::vector<Boost_Float> &normalization,
                      const size_t k, const mp_bitcnt_t precision)
  {
    ASSERT_EQUAL(z.size(), normalization.size());
    Scoped_Boost_Float_Precision scoped_precision(
      working_precision(z, normalization, k, precision));
    Boost_Float nz = -1;
    for(size_t i = 0; i < z.size(); ++i)
      nz += normalization[i] * Boost_Float(z[i]);
    return abs(nz);
  }
}

// Write z.txt, in the format of write_matrix(): "<N> 1", then one entry a line.
// NB: the choice of k is the inverse of the conversion
// Polynomial_Matrix_Program.objective -> Output_SDP.dual_objective_b
// in the Output_SDP constructor.
inline void write_z(const El::Matrix<El::BigFloat> &y,
                    const std::vector<Boost_Float> &normalization,
                    const std::filesystem::path &path)
{
  ASSERT_EQUAL(y.Width(), 1);
  const auto precision = El::gmp::Precision();
  const auto k = max_normalization_index(to_BigFloat_Vector(normalization));
  const auto z = z_normalization::z_entries(
    z_normalization::printed_entries(y), normalization, k, precision);

  std::ofstream os(path);
  os << z.size() << " 1\n";
  for(const auto &entry : z)
    os << entry << "\n";
  ASSERT(os.good(), "Error when writing to: ", path);

  // The writer gives |n.z - 1| ~ 2^-(P + 60): warn above 2^-P.
  // A warning only: z.txt is written, and the run goes on.
  try
    {
      const auto error
        = z_normalization::normalization_error(z, normalization, k, precision);
      if(error > ldexp(Boost_Float(1), -static_cast<int>(precision)))
        PRINT_WARNING(path, ": |n.z - 1| = ",
                      error.str(3, std::ios_base::scientific),
                      " > 2^-precision. z is not normalized as written.");
    }
  catch(std::exception &e)
    {
      PRINT_WARNING("Cannot check n.z == 1 for ", path, ": ", e.what());
    }
}
