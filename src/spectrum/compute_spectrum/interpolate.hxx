#pragma once

#include "pmp/Polynomial.hxx"
#include "sdpb_util/Boost_Float.hxx"

#include <El.hpp>

#include <vector>

// Lagrange basis for interpolation polynomials:
// l_i(x) = product_{j!=i}((x - x_j) / (x_i - x_j))
inline std::vector<Boost_Polynomial>
get_lagrange_basis(const std::vector<El::BigFloat> &sample_points)
{
  const auto num_points = sample_points.size();
  ASSERT(num_points > 0);
  std::vector<Boost_Polynomial> result;
  result.reserve(num_points);

  const auto x = to_Boost_Float_Vector(sample_points);

  const Boost_Float one = 1;
  for(size_t i = 0; i < num_points; ++i)
    {
      auto &poly = result.emplace_back(one);
      // poly = Product[(x - x_j) / (x_i - x_j), j != i]
      for(size_t j = 0; j < num_points; ++j)
        {
          if(j == i)
            continue;
          poly /= x.at(i) - x.at(j);
          poly *= Boost_Polynomial({-x.at(j), one});
        }
    }

  return result;
}

// Build interpolating polynomial through data points (x_i, y_i)
// using Lagrange formula:
// p(x) = sum_i(y_i * l_i(x))
// where l_i(x) = product_{j!=i}((x - x_j) / (x_i - x_j))
//
// Since we reuse the same x_values (aka sample_points) many times,
// we precompute Lagrange basis get_lagrange_basis(x_values)
// and pass it as a first argument instead of x_values themselves
inline Boost_Polynomial
interpolate(const std::vector<Boost_Polynomial> &lagrange_basis,
            const std::vector<El::BigFloat> &y_values)
{
  const auto num_points = y_values.size();
  ASSERT_EQUAL(num_points, lagrange_basis.size());

  Boost_Polynomial result;
  for(size_t i = 0; i < num_points; ++i)
    result += lagrange_basis.at(i) * to_Boost_Float(y_values.at(i));
  return result;
}

inline Boost_Polynomial interpolate(const std::vector<El::BigFloat> &x_values,
                                    const std::vector<El::BigFloat> &y_values)
{
  return interpolate(get_lagrange_basis(x_values), y_values);
}