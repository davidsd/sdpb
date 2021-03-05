#include "poles_prefactor.hxx"
#include "power_prefactor.hxx"
#include "Function.hxx"
#include "../sdp_read.hxx"

void convert_matrices_to_functions(
  const El::BigFloat &max_delta,
  const std::vector<Positive_Matrix_With_Prefactor> &matrices,
  std::vector<std::vector<std::vector<std::vector<Function>>>> &functions)
{
  const Boost_Float max_delta_mpfr([&]() {
    std::stringstream ss;
    set_stream_precision(ss);
    ss << max_delta;
    return ss.str();
  }());
  const Boost_Float pi(boost::math::constants::pi<Boost_Float>());

  for(auto &matrix_block : matrices)
    {
      const int64_t max_degree([&]() {
        int64_t result(0);
        for(auto &row : matrix_block.polynomials)
          for(auto &column : row)
            for(auto &poly : column)
              {
                result = std::max(result, poly.degree());
              }
        return result;
      }());

      const size_t N(max_degree+1);
      std::vector<El::BigFloat> x(N);
      for(size_t k(0); k < N; ++k)
        {
          Boost_Float x_unscaled(
                                 cos(pi * (2 * (N - 1 - k) + 1) / (2 * N)));
          x[k]=(El::BigFloat(to_string(x_unscaled)) + 1) * max_delta /2;
        }

      std::vector<std::vector<El::BigFloat>> weights(N);
      for(size_t n(0); n < N; ++n)
        {
          weights[n].reserve(N);
          for(size_t k(0); k < N; ++k)
            {
              Boost_Float weight_mpfr(
                                     2 * cos((n * pi * (2 * (N - 1 - k) + 1)) / (2 * N))
                                     / N);
              weights[n].emplace_back(to_string(weight_mpfr));
            }
        }
      
      std::vector<El::BigFloat> values(N);
      functions.emplace_back();
      auto &function_block(functions.back());
      for(auto &matrix_row : matrix_block.polynomials)
        {
          function_block.emplace_back();
          auto &function_row(function_block.back());
          for(auto &matrix_column : matrix_row)
            {
              function_row.emplace_back();
              auto &function_column(function_row.back());
              for(auto &poly : matrix_column)
                {
                  function_column.emplace_back();
                  auto &function(function_column.back());
                  function.max_delta = max_delta;
                  if(poly.degree() < max_degree)
                    {
                      function.infinity_value = 0;
                    }
                  else
                    {
                      function.infinity_value
                        = poly.coefficients.at(max_degree);
                    }
                  for(size_t k(0); k < N; ++k)
                    {
                      values[k] = poly(x[k]);
                      // We do not multiply by the Damped Rational
                      // prefactor.  It has poles, many near zero,
                      // which make any kind of approximation bad.
                      // Since positivity is not affected by the
                      // overall prefactor, we can omit it.
                    }

                  function.chebyshev_coeffs.reserve(N);
                  for(size_t n(0); n < N; ++n)
                    {
                      El::BigFloat coeff(0);
                      for(size_t k(0); k < N; ++k)
                        {
                          coeff += weights[n][k] * values[k];
                        }
                      function.chebyshev_coeffs.emplace_back(coeff);
                    }
                }
            }
        }
    }
}
