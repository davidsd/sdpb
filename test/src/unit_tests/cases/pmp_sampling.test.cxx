#include "pmp/Damped_Rational.hxx"
#include "sdpb_util/to_matrix.hxx"

#include <catch2/catch_amalgamated.hpp>

#include "unit_tests/util/util.hxx"

using Test_Util::REQUIRE_Equal::diff;

std::vector<Boost_Float>
sample_points(const size_t &num_points, const Damped_Rational &prefactor);

std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational);

std::array<El::Matrix<El::BigFloat>, 2>
sample_bilinear_basis(const std::vector<El::BigFloat> &sample_points,
                      const std::vector<El::BigFloat> &sample_scalings);

namespace
{
  Boost_Float exp_minus_one()
  {
    return
#if(BOOST_VERSION >= 106900)
      boost::math::constants::exp_minus_one<Boost_Float>();
#else
      Boost_Float(1) / boost::math::constants::e<Boost_Float>();
#endif
  }

  std::vector<El::BigFloat>
  to_BigFloat_Vector(const std::vector<Boost_Float> &input)
  {
    std::vector<El::BigFloat> output;
    output.reserve(input.size());
    for(const auto &x : input)
      {
        output.push_back(to_BigFloat(x));
      }

    return output;
  }

  void do_test_no_orig(const Damped_Rational &prefactor, const size_t degree)
  {
    INFO(
      "Trying to calculate sample points, scalings and bilinear bases to "
      "check if the code does not crash. Do not compare to expected values");

    CAPTURE(prefactor);

    const size_t num_points = degree + 1;
    CAPTURE(degree);
    CAPTURE(num_points);

    const auto boost_points = sample_points(num_points, prefactor);
    const auto points = to_BigFloat_Vector(boost_points);
    CAPTURE(points);

    const auto boost_scalings = sample_scalings(boost_points, prefactor);
    const auto scalings = to_BigFloat_Vector(boost_scalings);
    CAPTURE(scalings);

    const auto bilinear_bases = sample_bilinear_basis(points, scalings);
    CAPTURE(bilinear_bases[0]);
    CAPTURE(bilinear_bases[1]);
  }

  void do_test(
    const Damped_Rational &prefactor, const size_t degree,
    const std::vector<El::BigFloat> &expected_points,
    const std::vector<El::BigFloat> &expected_scalings,
    const std::array<El::Matrix<El::BigFloat>, 2> &expected_bilinear_bases,
    const int diff_precision)
  {
    CAPTURE(prefactor);

    const size_t num_points = degree + 1;
    CAPTURE(degree);
    CAPTURE(num_points);

    const auto boost_points = sample_points(num_points, prefactor);
    const auto points = to_BigFloat_Vector(boost_points);
    CAPTURE(points);

    const auto boost_scalings = sample_scalings(boost_points, prefactor);
    const auto scalings = to_BigFloat_Vector(boost_scalings);
    CAPTURE(scalings);

    const auto bilinear_bases = sample_bilinear_basis(points, scalings);
    CAPTURE(bilinear_bases[0]);
    CAPTURE(bilinear_bases[1]);

    DIFF_PREC(points, expected_points, diff_precision);
    DIFF_PREC(scalings, expected_scalings, diff_precision);
    DIFF_PREC(bilinear_bases[0], expected_bilinear_bases[0], diff_precision);
    DIFF_PREC(bilinear_bases[1], expected_bilinear_bases[1], diff_precision);
  }
}

TEST_CASE("pmp_sampling")
{
  if(El::mpi::Rank() != 0)
    return;

  const int diff_precision = 16;

  SECTION("Compare with Mathematica")
  {
    SECTION("exp(-x), degree=4")
    {
      const Damped_Rational prefactor{1, exp_minus_one(), {}};
      const size_t degree = 4;

      // The numbers below are calculated in Mathematica.
      // Truncated to 5 decimal digits, i.e. precision = 16.

      const std::vector<El::BigFloat> points{0.061812, 0.56588, 1.6319, 3.4239,
                                             6.4864};
      const std::vector<El::BigFloat> scalings{0.94006, 0.56786, 0.19555,
                                               0.032586, 0.0015240};
      std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases;

      // NB: some rows of bilinear_bases in C++ implementation
      // have different sign compared to Mathematica.
      // This doesn't really matter:
      // one can choose the signs arbitrarily when performing SVD decomposition of a matrix.
      // To compensate for this, we multiply rows by -1 when necessary to pass the test.

      bilinear_bases[0] = to_matrix<El::BigFloat>(
        {{0.059050, 0.15519, 0.44912, 0.70680, 0.52072},
         {-0.72208, -0.58512, -0.23133, 0.15942, 0.23941},
         {0.58115, -0.27484, -0.59714, 0.038956, 0.47816}});
      bilinear_bases[0](2, El::ALL) *= -1;

      bilinear_bases[1] = to_matrix<El::BigFloat>(
        {{0.065390, 0.29723, 0.59855, 0.65441, 0.34766},
         {0.40713, 0.71945, 0.21543, -0.37109, -0.36404}});
      bilinear_bases[1](0, El::ALL) *= -1;

      do_test(prefactor, degree, points, scalings, bilinear_bases, diff_precision);
    }
  }

  SECTION("Crash tests")
  {
    for(const std::vector<Boost_Float> &poles :
        std::vector<std::vector<Boost_Float>>{
          {}, {0}, {0, 0}, {-1}, {0, -1}, {0, -1, -2}})
      {
        DYNAMIC_SECTION("poles=" << poles)
        {
          const size_t degree = GENERATE(0, 1, 2, 10);
          DYNAMIC_SECTION("degree=" << degree)
          {
            const Damped_Rational prefactor{1, exp_minus_one(), poles};
            do_test_no_orig(prefactor, degree);
          }
        }
      }
  }
}
