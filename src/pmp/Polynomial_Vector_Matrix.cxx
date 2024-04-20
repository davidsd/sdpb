#include "Polynomial_Vector_Matrix.hxx"

#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"

#include <boost/math/constants/constants.hpp>

std::vector<Boost_Float>
sample_points(const size_t &num_points, const Damped_Rational &prefactor);

std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational);

std::array<Polynomial_Vector, 2>
bilinear_basis(const std::vector<El::BigFloat> &sample_points,
               const std::vector<El::BigFloat> &sample_scalings);

namespace
{
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

  std::vector<Boost_Float>
  to_Boost_Float_Vector(const std::vector<El::BigFloat> &input)
  {
    std::vector<Boost_Float> output;
    output.reserve(input.size());
    for(const auto &x : input)
      {
        output.push_back(to_Boost_Float(x));
      }

    return output;
  }
}

namespace
{
  int64_t get_max_degree(const El::Matrix<Polynomial_Vector> &pvm)
  {
    int64_t max_degree = 0;
    for(int i = 0; i < pvm.Height(); ++i)
      for(int j = 0; j < pvm.Width(); ++j)
        for(auto &polynomial : pvm(i, j))
          max_degree = std::max(max_degree, polynomial.degree());
    return max_degree;
  }

  // We use XXX_or_default() functions instead of std::optional::value_or()
  // to prevent unnecessary computation of default value.

  Damped_Rational
  prefactor_or_default(const std::optional<Damped_Rational> &prefactor_opt)
  {
    if(prefactor_opt.has_value())
      return prefactor_opt.value();

    // NB: we cannot initialize exp_minus_one as a global constant,
    // because we don't know required MPFR precision at that moment.
    // Thus, we do it inside the function.
    const Boost_Float exp_minus_one =
#if(BOOST_VERSION >= 106900)
      boost::math::constants::exp_minus_one<Boost_Float>();
#else
      Boost_Float(1) / boost::math::constants::e<Boost_Float>();
#endif

    // Default prefactor is e^-x
    const Damped_Rational exp_minus_x = {1, exp_minus_one, {}};
    return exp_minus_x;
  }

  std::vector<El::BigFloat> sample_points_or_default(
    const std::optional<std::vector<El::BigFloat>> &sample_points_opt,
    const Damped_Rational &prefactor, const size_t max_degree)
  {
    if(sample_points_opt.has_value())
      return sample_points_opt.value();

    return to_BigFloat_Vector(sample_points(max_degree + 1, prefactor));
  }

  std::vector<El::BigFloat> sample_scalings_or_default(
    const std::optional<std::vector<El::BigFloat>> &sample_scalings_opt,
    const std::vector<El::BigFloat> &sample_points,
    const Damped_Rational &damped_rational)
  {
    if(sample_scalings_opt.has_value())
      return sample_scalings_opt.value();
    return to_BigFloat_Vector(
      sample_scalings(to_Boost_Float_Vector(sample_points), damped_rational));
  }

  std::array<Polynomial_Vector, 2> bilinear_basis_or_default(
    const std::optional<Polynomial_Vector> &bilinear_basis_opt,
    const std::vector<El::BigFloat> &sample_points,
    const std::vector<El::BigFloat> &sample_scalings)
  {
    if(bilinear_basis_opt.has_value())
      {
        const auto &basis_all = bilinear_basis_opt.value();

        std::array<Polynomial_Vector, 2> basis;
        const size_t degree = sample_points.size() - 1;
        const size_t delta1 = degree / 2;

        // Check input size.
        // delta1 >= delta2 always, so we need to compare only with (delta1 + 1)
        if(basis_all.size() < delta1 + 1)
          {
            RUNTIME_ERROR("PMP: bilinearBases size=", basis_all.size(),
                          ", required at least ", delta1 + 1);
          }
        if(basis_all.size() > delta1 + 1)
          {
            PRINT_WARNING("PMP: bilinearBases size=", basis_all.size(),
                          " is too large, only the first ", delta1 + 1,
                          " polynomials will be used");
          }

        basis[0] = Polynomial_Vector(basis_all.begin(),
                                     basis_all.begin() + delta1 + 1);
        if(degree > 0)
          {
            const size_t delta2 = (degree + 1) / 2 - 1;
            basis[1] = Polynomial_Vector(basis_all.begin(),
                                         basis_all.begin() + delta2 + 1);
            // otherwise keep basis[1] empty
          }
        return basis;
      }

    return bilinear_basis(sample_points, sample_scalings);
  }
}

Polynomial_Vector_Matrix::Polynomial_Vector_Matrix(
  const El::Matrix<Polynomial_Vector> &polynomials,
  const std::optional<Damped_Rational> &prefactor_opt,
  const std::optional<Damped_Rational> &reduced_prefactor_opt,
  const std::optional<std::vector<El::BigFloat>> &sample_points_opt,
  const std::optional<std::vector<El::BigFloat>> &sample_scalings_opt,
  const std::optional<std::vector<El::BigFloat>> &reduced_sample_scalings_opt,
  const std::optional<Polynomial_Vector> &bilinear_basis_opt)
{
  this->polynomials = polynomials;
  const auto prefactor = prefactor_or_default(prefactor_opt);

  const auto reduced_prefactor = [&] {
    if(reduced_prefactor_opt.has_value())
      {
        ASSERT(prefactor_opt.has_value());
        return prefactor_or_default(reduced_prefactor_opt);
      }
    return prefactor;
  }();

  if(reduced_prefactor.poles.size() > prefactor.poles.size())
    {
      PRINT_WARNING("reducedPrefactor has more poles than prefactor, the "
                    "number of sample points will be increased!\n\t",
                    DEBUG_STRING(prefactor.poles.size()),
                    DEBUG_STRING(reduced_prefactor.poles.size()));
    }

  const int64_t max_degree = get_max_degree(polynomials)
                             + reduced_prefactor.poles.size()
                             - prefactor.poles.size();
  ASSERT(max_degree >= 0, DEBUG_STRING(max_degree), DEBUG_STRING(prefactor),
         DEBUG_STRING(reduced_prefactor));

  sample_points = sample_points_or_default(sample_points_opt,
                                           reduced_prefactor, max_degree);

  sample_scalings = sample_scalings_or_default(sample_scalings_opt,
                                               this->sample_points, prefactor);
  reduced_sample_scalings = [&] {
    if(reduced_sample_scalings_opt.has_value()
       || reduced_prefactor_opt.has_value())
      {
        return sample_scalings_or_default(
          reduced_sample_scalings_opt, this->sample_points, reduced_prefactor);
      }
    return sample_scalings;
  }();

  bilinear_basis = bilinear_basis_or_default(bilinear_basis_opt, sample_points,
                                             reduced_sample_scalings);

  validate(max_degree);
}

void Polynomial_Vector_Matrix::validate(const int64_t max_degree) const
{
  ASSERT_EQUAL(sample_points.size(), max_degree + 1);
  ASSERT_EQUAL(reduced_sample_scalings.size(), sample_points.size());

  const size_t delta1 = max_degree / 2;
  ASSERT_EQUAL(bilinear_basis[0].size(), delta1 + 1, DEBUG_STRING(max_degree));
  if(max_degree == 0)
    {
      ASSERT_EQUAL(bilinear_basis[1].size(), 0);
    }
  else
    {
      const size_t delta2 = (max_degree + 1) / 2 - 1;
      ASSERT_EQUAL(bilinear_basis[1].size(), delta2 + 1,
                   DEBUG_STRING(max_degree));
    }

  ASSERT_EQUAL(polynomials.Height(), polynomials.Width());

  // TODO check if this is fast enough
  for(int i = 0; i < polynomials.Height(); ++i)
    for(int j = 0; j < polynomials.Width(); ++j)
      {
        if(i == j)
          continue;
        ASSERT_EQUAL(polynomials(i, j), polynomials(j, i), DEBUG_STRING(i),
                     DEBUG_STRING(j));
      }
}
