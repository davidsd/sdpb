#include "../Function_State.hxx"
#include "sdpb_util/Boost_Float.hxx"

#include <boost/math/constants/constants.hpp>

void Function_State::json_end_object()
{
  if(parsing_max_delta)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + max_delta_state.name + "'.");
    }
  else if(parsing_epsilon_value)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + epsilon_value_state.name + "'.");
    }
  else if(parsing_infinity_value)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + infinity_value_state.name + "'.");
    }
  else if(parsing_chebyshev_values)
    {
      throw std::runtime_error(
        "Invalid input file.  Unexpected object end inside '" + name + "."
        + chebyshev_values_state.name + "'.");
    }

  // Convert from sampled values to chebyshev coefficients.
  const Boost_Float pi(boost::math::constants::pi<Boost_Float>());
  auto &values(chebyshev_values_state.value);
  const size_t N(values.size());
  value.chebyshev_coeffs.resize(0);
  value.chebyshev_coeffs.reserve(N);
  std::vector<Boost_Float> coeffs(N, Boost_Float(0));
  for(size_t n(0); n < N; ++n)
    {
      Boost_Float coeff(0);
      for(size_t k(0); k < N; ++k)
        {
          coeff += 2 * cos((n * pi * (2 * (N - 1 - k) + 1)) / (2 * N))
                   * values[k] / N;
        }
      value.chebyshev_coeffs.emplace_back(to_string(coeff));
    }
  inside = false;
}
