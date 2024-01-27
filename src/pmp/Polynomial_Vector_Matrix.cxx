#include "Polynomial_Vector_Matrix.hxx"

#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"

#include <boost/math/constants/constants.hpp>

std::vector<Boost_Float> sample_points(const size_t &num_points);

std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational);

Polynomial_Vector bilinear_basis(const Damped_Rational &damped_rational,
                                 const size_t &half_max_degree);

namespace
{
  std::vector<El::BigFloat>
  to_BigFloat_Vector(const std::vector<Boost_Float> &input)
  {
    std::vector<El::BigFloat> output;
    output.reserve(input.size());

    std::stringstream ss;
    set_stream_precision(ss);

    for(const auto &x : input)
      {
        ss.str("");
        ss << x;
        output.emplace_back(ss.str());
      }

    return output;
  }

  std::vector<Boost_Float>
  to_Boost_Float_Vector(const std::vector<El::BigFloat> &input)
  {
    std::vector<Boost_Float> output;
    output.reserve(input.size());

    std::stringstream ss;
    set_stream_precision(ss);

    for(const auto &x : input)
      {
        ss.str("");
        ss << x;
        output.emplace_back(ss.str());
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

  // Default prefactor is e^-x
  const Damped_Rational exp_minus_x
    = {1, boost::math::constants::exp_minus_one<Boost_Float>(), {}};

  Damped_Rational
  prefactor_or_default(const std::optional<Damped_Rational> &prefactor_opt)
  {
    return prefactor_opt.has_value() ? prefactor_opt.value() : exp_minus_x;
  }

  std::vector<El::BigFloat> sample_points_or_default(
    const std::optional<std::vector<El::BigFloat>> &sample_points_opt,
    const size_t max_degree)
  {
    if(sample_points_opt.has_value())
      return sample_points_opt.value();

    return to_BigFloat_Vector(sample_points(max_degree + 1));
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
  Polynomial_Vector bilinear_basis_or_default(
    const std::optional<Polynomial_Vector> &bilinear_basis_opt,
    const Damped_Rational &damped_rational, int64_t max_degree)
  {
    if(bilinear_basis_opt.has_value())
      return bilinear_basis_opt.value();

    return bilinear_basis(damped_rational, max_degree / 2);
  }
}

Polynomial_Vector_Matrix::Polynomial_Vector_Matrix(
  const El::Matrix<Polynomial_Vector> &polynomials,
  const std::optional<Damped_Rational> &prefactor_opt,
  const std::optional<std::vector<El::BigFloat>> &sample_points_opt,
  const std::optional<std::vector<El::BigFloat>> &sample_scalings_opt,
  const std::optional<Polynomial_Vector> &bilinear_basis_opt)
{
  const auto max_degree = get_max_degree(polynomials);

  this->polynomials = polynomials;
  const auto prefactor = prefactor_or_default(prefactor_opt);
  sample_points = sample_points_or_default(sample_points_opt, max_degree);
  sample_scalings = sample_scalings_or_default(sample_scalings_opt,
                                               this->sample_points, prefactor);
  bilinear_basis
    = bilinear_basis_or_default(bilinear_basis_opt, prefactor, max_degree);

  validate(max_degree);
}

void Polynomial_Vector_Matrix::validate(const int64_t max_degree) const
{
  ASSERT(sample_points.size() == max_degree + 1,
         "sample_points.size()=", sample_points.size(), ", expected ",
         max_degree + 1);
  ASSERT(sample_scalings.size() == sample_points.size(),
         "sample_scalings.size()=", sample_scalings.size(), ", expected ",
         sample_points.size());
  ASSERT(bilinear_basis.size() == max_degree / 2 + 1,
         "bilinear_basis.size()=", bilinear_basis.size(), ", expected ",
         max_degree / 2 + 1);

  ASSERT(polynomials.Height() == polynomials.Width(),
         "Height()=", polynomials.Height(), " Width()=", polynomials.Width());

  // TODO check if this is fast enough
  for(int i = 0; i < polynomials.Height(); ++i)
    for(int j = 0; j < polynomials.Width(); ++j)
      {
        if(i == j)
          continue;
        ASSERT(polynomials(i, j) == polynomials(j, i), "i=", i, " j=", j);
      }
}
