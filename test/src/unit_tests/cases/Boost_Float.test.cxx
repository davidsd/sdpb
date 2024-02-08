#include <catch2/catch_amalgamated.hpp>

#include "test_util/diff.hxx"
#include "sdpb_util/Boost_Float.hxx"

using Test_Util::REQUIRE_Equal::diff;

TEST_CASE("mpfr2mpf")
{
  Boost_Float::default_precision(3);
  for(const El::BigFloat &input :
      std::initializer_list<El::BigFloat>{1.0, -123.45, 56434.4536645})
    {
      CAPTURE(input);
      auto output = to_BigFloat(to_Boost_Float(input));
      DIFF(input, output);
    }

  El::BigFloat a = 13423.1324, b = -142343e120;
  Boost_Float a_2 = to_Boost_Float(a), b_2 = to_Boost_Float(b);
  Boost_Float::default_precision(3);
  DIFF(a * b, to_BigFloat(a_2 * b_2));
}
