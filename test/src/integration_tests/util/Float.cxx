#include "Float.hxx"

#include "diff.hxx"

#include <El.hpp>

namespace
{
  // Copied from El::gmp::SetPrecision()
  // We remove MPI code and
  // allow for changing precision many times (at our own risk)
  // - this is needed for our tests.
  void El_gmp_set_precision(mp_bitcnt_t prec)
  {
    El::gmp::num_limbs
      = (std::max(prec, static_cast<mp_bitcnt_t>(53)) + 2 * GMP_NUMB_BITS - 1)
          / GMP_NUMB_BITS
        + 1;
    mpf_set_default_prec(prec);
  }
}

Float_Binary_Precision::Float_Binary_Precision(unsigned int binary_precision,
                                               unsigned int diff_precision)
    : old_binary_precision(binary_precision),
      old_diff_precision(diff_precision)
{
  El_gmp_set_precision(binary_precision);
  Test_Util::REQUIRE_Equal::diff_precision = diff_precision;
}
Float_Binary_Precision::~Float_Binary_Precision()
{
  El_gmp_set_precision(old_binary_precision);
  Test_Util::REQUIRE_Equal::diff_precision = old_diff_precision;
}
