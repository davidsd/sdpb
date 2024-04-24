#include "pmp/Damped_Rational.hxx"
#include "pmp/Polynomial.hxx"
#include "sdpb_util/to_matrix.hxx"

#include <catch2/catch_amalgamated.hpp>

#include "unit_tests/util/util.hxx"

using Test_Util::REQUIRE_Equal::diff;

std::vector<Boost_Float>
sample_points(const size_t &num_points, const Damped_Rational &prefactor);

std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational);

std::array<Polynomial_Vector, 2>
bilinear_basis(const std::vector<El::BigFloat> &sample_points,
               const std::vector<El::BigFloat> &sample_scalings);

std::array<El::Matrix<El::BigFloat>, 2>
sample_bilinear_basis(const std::array<Polynomial_Vector, 2> &bilinear_basis,
                      const std::vector<El::BigFloat> &sample_points,
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

    const auto bilinear_basis_arr = bilinear_basis(points, scalings);
    std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases
      = sample_bilinear_basis(bilinear_basis_arr, points, scalings);
    CAPTURE(bilinear_basis_arr[0]);
    CAPTURE(bilinear_basis_arr[1]);
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

    const auto bilinear_basis_arr = bilinear_basis(points, scalings);
    std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases
      = sample_bilinear_basis(bilinear_basis_arr, points, scalings);
    CAPTURE(bilinear_basis_arr[0]);
    CAPTURE(bilinear_basis_arr[1]);

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
      // TODO: Mathematica's bilinear bases are computed via SVD,
      // we should update them to Cholesky version

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
        {{0.735538, 0.571674, 0.335472, 0.136945, 0.0296154},
         {-0.454533, 0.0809192, 0.586352, 0.609108, 0.268386},
         {0.339114, -0.329580, -0.394268, 0.369777, 0.695842}});

      bilinear_bases[1] = to_matrix<El::BigFloat>(
        {{0.266195, 0.625991, 0.623829, 0.368860, 0.109794},
         {-0.314913, -0.462694, 0.124527, 0.655667, 0.491262}});

      do_test(prefactor, degree, points, scalings, bilinear_bases,
              diff_precision);
    }

    SECTION("exp(-x)/x/(x+1), degree=4")
    {
      const Damped_Rational prefactor{1, exp_minus_one(), {0, -1}};
      const size_t degree = 4;

      // The numbers below are calculated in Mathematica.
      // Truncated to 5 decimal digits, i.e. precision = 16.
      // TODO: Mathematica's bilinear bases are computed via SVD,
      // we should update them to Cholesky version

      const std::vector<El::BigFloat> points{0, 0.0501905, 0.490170, 1.56871,
                                             3.82960};
      const std::vector<El::BigFloat> scalings{1.0000e16, 18.0432, 0.838570,
                                               0.0516961, 0.00117426};
      std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases;

      // NB: some rows of bilinear_bases in C++ implementation
      // have different sign compared to Mathematica.
      // This doesn't really matter:
      // one can choose the signs arbitrarily when performing SVD decomposition of a matrix.
      // To compensate for this, we multiply rows by -1 when necessary to pass the test.

      bilinear_bases[0] = to_matrix<El::BigFloat>(
        {{1.00000, 4.247724e-8, 9.157346e-9, 2.273678e-9, 3.426744e-10},
         {-2.241430e-8, 0.3407876, 0.7175002, 0.5701354, 0.2097686},
         {1.771584e-8, -0.3631303, -0.3850530, 0.4332262, 0.7295105}});

      bilinear_bases[1] = to_matrix<El::BigFloat>(
        {{0, 0.8036324, 0.5414187, 0.2404866, 0.05663024},
         {0, -0.4294390, 0.2667574, 0.7239642, 0.4693597}});

      do_test(prefactor, degree, points, scalings, bilinear_bases,
              diff_precision);
    }

    SECTION("prefactor=1, degree=0")
    {
      INFO("Special case: constant prefactor and constant constraints. "
           "In that case, we can choose any sample point. "
           "SDPB chooses 0.");
      const size_t degree = 0;
      const Damped_Rational prefactor{1, 1, {}};

      const std::vector<El::BigFloat> points{0};
      const std::vector<El::BigFloat> scalings{1};
      std::array<El::Matrix<El::BigFloat>, 2> bilinear_bases;
      bilinear_bases[0] = to_matrix<El::BigFloat>({{1}});
      bilinear_bases[1] = El::Matrix<El::BigFloat>(0, 1);

      do_test(prefactor, degree, points, scalings, bilinear_bases,
              diff_precision);
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
