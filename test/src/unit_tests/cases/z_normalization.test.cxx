#include <catch2/catch_amalgamated.hpp>
#include "sdpb/z_normalization.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "test_util/diff.hxx"

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// z.txt's eliminated component z_k is formed so that the written decimals
// satisfy n.z == 1 (sdpb/z_normalization.hxx).

using Test_Util::REQUIRE_Equal::diff;

namespace
{
  // n as SDP::normalization reads it: all digits of the decimals
  std::vector<Boost_Float> values(const std::vector<std::string> &decimals)
  {
    Scoped_Boost_Float_Precision precision(1000);
    std::vector<Boost_Float> result;
    for(const auto &s : decimals)
      result.emplace_back(s);
    return result;
  }

  // log2|n.z - 1| for the values of n (as SDPB holds them) and the decimals
  // of z, at 5000 bits (exact enough here)
  int log2_error(const std::vector<Boost_Float> &n,
                 const std::vector<std::string> &z)
  {
    REQUIRE(n.size() == z.size());
    Scoped_Boost_Float_Precision precision(5000);
    Boost_Float nz = -1;
    for(size_t i = 0; i < n.size(); ++i)
      nz += n[i] * Boost_Float(z[i]);
    if(nz == 0)
      return std::numeric_limits<int>::min();
    int exponent;
    frexp(nz, &exponent);
    return exponent;
  }

  const std::string sqrt_half
    = "0.7071067811865475244008443621048490392848359376884740365883398689953"
      "6623923105351942519376716382078636750692311545614851246241802792536860"
      "632206074854996791570661133296375279637789997525057639103028574";
  const std::string inv_pi
    = "0.3183098861837906715377675267450287240689192914809128974953346881177"
      "9359526845307018022760553250617191214568545351591607378582369222915730"
      "575593482146339967845847993387481815514615549279385061537743478";
  const std::string third = "-0.333333333333333333333333333333333333333333333"
                            "33333333333333333333333333333333333333333";
}

TEST_CASE("z_normalization")
{
  SECTION("huge y at low precision: written z has n.z == 1")
  {
    // y as SDPB at 64 bits would print it (~20 digits); |n_i y_i| ~ 1e84
    const std::vector<std::string> y{"-9.5179409895920326840e+83",
                                     "4.4741791903159498881e+84"};
    const std::vector<std::string> n{third, sqrt_half, inv_pi};
    for(const mp_bitcnt_t precision : {64, 128, 832})
      for(const size_t k : {0, 1, 2})
        {
          DYNAMIC_SECTION("precision=" << precision << " k=" << k)
          {
            const auto z
              = z_normalization::z_entries(y, values(n), k, precision);
            REQUIRE(z.size() == 3);
            // the other components are y's strings
            auto z_without_k = z;
            z_without_k.erase(z_without_k.begin() + k);
            DIFF(z_without_k, y);
            const int error = log2_error(values(n), z);
            CAPTURE(z[k], error);
            REQUIRE(error <= -static_cast<int>(precision) - 32);
            // the self-check agrees
            REQUIRE(z_normalization::normalization_error(z, values(n), k,
                                                         precision)
                    <= ldexp(Boost_Float(1), -static_cast<int>(precision)));
          }
        }
  }

  SECTION("z_k formed at the solver's precision fails the same check")
  {
    // The old way: z_0 = (1 - n_1 y_1 - n_2 y_2) / n_0 at 64 bits,
    // printed with 20 digits. Shows that the check above can fail.
    const std::vector<std::string> y{"-9.5179409895920326840e+83",
                                     "4.4741791903159498881e+84"};
    const std::vector<std::string> n{sqrt_half, third, inv_pi};
    Boost_Float z0;
    {
      Scoped_Boost_Float_Precision precision(64);
      z0 = (1 - Boost_Float(n[1]) * Boost_Float(y[0])
            - Boost_Float(n[2]) * Boost_Float(y[1]))
           / Boost_Float(n[0]);
    }
    std::ostringstream os;
    os.precision(20);
    os << z0;
    const std::vector<std::string> z{os.str(), y[0], y[1]};
    REQUIRE(log2_error(values(n), z) > 100);
    // and the self-check flags it
    REQUIRE(z_normalization::normalization_error(z, values(n), 0, 64)
            > ldexp(Boost_Float(1), -64));
  }

  SECTION("one component (N = 1)")
  {
    const std::vector<std::string> n{inv_pi};
    const auto z = z_normalization::z_entries({}, values(n), 0, 64);
    REQUIRE(z.size() == 1);
    REQUIRE(log2_error(values(n), z) <= -64 - 32);
  }
}
