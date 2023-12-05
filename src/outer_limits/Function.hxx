#pragma once

#include <El.hpp>

#include "sdpb_util/ostream/ostream_vector.hxx"

struct Function
{
  El::BigFloat max_delta, epsilon_value, infinity_value;
  std::vector<El::BigFloat> chebyshev_coeffs;

  El::BigFloat eval(const El::BigFloat &epsilon,
                    const El::BigFloat &infinity,
                    const El::BigFloat &x) const;
};

inline std::ostream & operator<<(std::ostream &os, const Function &f)
{
  os << "{\n  \"max_delta\": "
     << f.max_delta
     << ",\n  \"epsilon_value\": "
     << f.epsilon_value
     << ",\n  \"infinity_value\": "
     << f.infinity_value
     << ",\n  \"chebyshev_coeffs\": "
     << f.chebyshev_coeffs
     << "\n}";
  return os;
}
  
