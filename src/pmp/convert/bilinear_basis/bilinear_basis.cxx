#include "pmp_read/pmp_read.hxx"

void precompute(
  const Boost_Float &base, std::vector<Boost_Float> &sorted_poles,
  std::vector<std::pair<std::vector<Boost_Float>::const_iterator,
                        std::vector<Boost_Float>::const_iterator>>
    &equal_ranges,
  std::vector<int64_t> &lengths, std::vector<Boost_Float> &products,
  std::vector<std::vector<Boost_Float>> &integral_matrix);

Boost_Float bilinear_form(
  const Damped_Rational &damped_rational,
  const std::vector<Boost_Float> &sorted_poles,
  const std::vector<std::pair<std::vector<Boost_Float>::const_iterator,
                              std::vector<Boost_Float>::const_iterator>>
    &equal_ranges,
  const std::vector<int64_t> &lengths,
  const std::vector<Boost_Float> &products,
  const std::vector<std::vector<Boost_Float>> &integral_matrix,
  const int64_t &m);

Polynomial_Vector bilinear_basis(const Damped_Rational &damped_rational,
                                 const size_t &half_max_degree)
{
  // Exit early if damped_rational is a constant
  if(damped_rational.is_constant())
    {
      Polynomial_Vector result;
      result.emplace_back(1, to_BigFloat(1 / sqrt(damped_rational.constant)));
      return result;
    }

  Polynomial polynomial(half_max_degree, 1);

  std::vector<Boost_Float> sorted_poles(damped_rational.poles);
  std::vector<std::pair<std::vector<Boost_Float>::const_iterator,
                        std::vector<Boost_Float>::const_iterator>>
    equal_ranges;
  std::vector<int64_t> lengths;
  std::vector<Boost_Float> products;
  std::vector<std::vector<Boost_Float>> integral_matrix;
  precompute(damped_rational.base, sorted_poles, equal_ranges, lengths,
             products, integral_matrix);

  std::vector<El::BigFloat> bilinear_table;
  for(int64_t m = 0; m <= int64_t(2 * half_max_degree); ++m)
    {
      bilinear_table.emplace_back(
        to_BigFloat(bilinear_form(damped_rational, sorted_poles, equal_ranges,
                                  lengths, products, integral_matrix, m)));
    }

  El::Matrix<El::BigFloat> anti_band_matrix(half_max_degree + 1,
                                            half_max_degree + 1);
  for(size_t i = 0; i < half_max_degree + 1; ++i)
    for(size_t j = i; j < half_max_degree + 1; ++j)
      {
        anti_band_matrix.Set(half_max_degree - j, j - i,
                             bilinear_table[half_max_degree - i]);
        if(i > 0)
          {
            anti_band_matrix.Set(half_max_degree - (j - i), j,
                                 bilinear_table[half_max_degree + i]);
          }
      }

  Cholesky(El::UpperOrLowerNS::UPPER, anti_band_matrix);

  El::Matrix<El::BigFloat> basis(half_max_degree + 1, half_max_degree + 1);

  El::Fill(basis, El::BigFloat(0));
  El::FillDiagonal(basis, El::BigFloat(1));

  El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::UPPER,
           El::OrientationNS::NORMAL, El::UnitOrNonUnit::NON_UNIT,
           El::BigFloat(1), anti_band_matrix, basis);

  Polynomial_Vector result(basis.Height());
  for(int64_t row = 0; row < basis.Height(); ++row)
    for(int64_t column = 0; column < basis.Width(); ++column)
      {
        if(basis(column, row) != El::BigFloat(0))
          {
            result[row].coefficients.resize(column + 1, 0);
            result[row].coefficients[column] = basis(column, row);
          }
      }
  return result;
}
