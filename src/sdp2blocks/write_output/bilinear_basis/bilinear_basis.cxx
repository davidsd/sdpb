#include "../../Damped_Rational.hxx"
#include "../../../Polynomial.hxx"
#include "../../../Timers.hxx"

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

std::vector<Polynomial>
bilinear_basis(const Damped_Rational &damped_rational,
               const size_t &half_max_degree, const std::string &timer_prefix,
               Timers &timers)
{
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
      const std::string bilinear_form_timer_name(
        timer_prefix + ".bilinear_form_" + std::to_string(m));
      auto &bilinear_form_timer(
        timers.add_and_start(bilinear_form_timer_name));

      // Using a stringstream seems to the best way to convert between
      // MPFR and GMP.  It may lose a bit or two since string
      // conversion is not sufficient for round-tripping.
      std::stringstream ss;
      ss.precision(Boost_Float::default_precision());
      ss << bilinear_form(damped_rational, sorted_poles, equal_ranges, lengths,
                          products, integral_matrix, m);
      bilinear_table.emplace_back(ss.str());
      bilinear_form_timer.stop();
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

  std::vector<Polynomial> result(basis.Height());
  for(int64_t row=0; row<basis.Height(); ++row)
    for(int64_t column=0; column<basis.Width(); ++column)
        {
          if(basis(column,row)!=El::BigFloat(0))
            {
              result[row].coefficients.resize(column+1,0);
              result[row].coefficients[column]=basis(column,row);
            }
      }
  return result;
}
