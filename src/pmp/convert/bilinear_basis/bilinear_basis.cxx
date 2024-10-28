#include "pmp_read/pmp_read.hxx"

namespace
{
  Polynomial_Vector compute_orthogonal_polynomials(
    const std::vector<El::BigFloat> &bilinear_table)
  {
    ASSERT_EQUAL(bilinear_table.size() % 2, 1,
                 DEBUG_STRING(bilinear_table.size()));

    const size_t delta = bilinear_table.size() / 2;

    El::Matrix<El::BigFloat> anti_band_matrix(delta + 1, delta + 1);
    for(size_t i = 0; i < delta + 1; ++i)
      for(size_t j = i; j < delta + 1; ++j)
        {
          anti_band_matrix.Set(delta - j, j - i, bilinear_table[delta - i]);
          if(i > 0)
            {
              anti_band_matrix.Set(delta - (j - i), j,
                                   bilinear_table[delta + i]);
            }
        }

    Cholesky(El::UpperOrLowerNS::UPPER, anti_band_matrix);
    // print warning if matrix is ill-conditioned
    {
      El::BigFloat min = anti_band_matrix.Get(0, 0);
      El::BigFloat max = min;
      for(int i = 0; i < anti_band_matrix.Height(); ++i)
        {
          const auto &val = anti_band_matrix.CRef(i, i);
          min = std::min(min, val);
          max = std::max(max, val);
        }
      // Condition number for Cholesky matrix can be estimated as a squared ratio of maximal and minimal elements on its diagonal,
      // see https://scicomp.stackexchange.com/questions/32762/cholmod-condition-number-estimate
      const auto condition_number = max * max / min / min;
      // 2^{precision/2}
      const auto threshold = El::BigFloat(1) <<= El::gmp::Precision() / 2;
      if(condition_number > threshold)
        {
          PRINT_WARNING("bilinear bases: anti_band_matrix is ill-conditioned, "
                        "this may reduce SDPB accuracy: ",
                        DEBUG_STRING(condition_number));
        }
    }

    El::Matrix<El::BigFloat> basis(delta + 1, delta + 1);

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
}

std::array<Polynomial_Vector, 2>
bilinear_basis(const std::vector<El::BigFloat> &sample_points,
               const std::vector<El::BigFloat> &sample_scalings)
{
  const size_t degree = sample_points.size() - 1;
  if(degree == 0)
    {
      Polynomial_Vector basis_0, basis_1;
      basis_0.push_back(Polynomial(1, 1));
      return {basis_0, basis_1};
    }

  std::vector<El::BigFloat> bilinear_table_all(degree + 1, 0);

  for(size_t i = 0; i < sample_points.size(); ++i)
    {
      El::BigFloat x_pow_n = 1;
      for(size_t n = 0; n < bilinear_table_all.size(); ++n)
        {
          bilinear_table_all.at(n) += x_pow_n * sample_scalings.at(i);
          x_pow_n *= sample_points.at(i);
        }
    }

  std::array<std::vector<El::BigFloat>, 2> bilinear_tables;
  const size_t delta1 = degree / 2;
  const size_t delta2 = (degree + 1) / 2 - 1;

  bilinear_tables[0] = std::vector(
    bilinear_table_all.begin(), bilinear_table_all.begin() + 2 * delta1 + 1);

  bilinear_tables[1]
    = std::vector(bilinear_table_all.begin() + 1,
                  bilinear_table_all.begin() + 1 + 2 * delta2 + 1);

  return {compute_orthogonal_polynomials(bilinear_tables[0]),
          compute_orthogonal_polynomials(bilinear_tables[1])};
}
