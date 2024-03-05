#include "sdpb_util/Boost_Float.hxx"

#include <El.hpp>

#include <boost/math/constants/constants.hpp>

#include <vector>

// Rescaled Laguerre
std::vector<Boost_Float> sample_points(const size_t &num_points)
{
  std::vector<Boost_Float> result;
  result.reserve(num_points);

  Boost_Float rho_crossing(3 - 2 * sqrt(Boost_Float(2))),
    constant(-boost::math::constants::pi_sqr<Boost_Float>()
                    / (64 * num_points * log(rho_crossing)));

  for(size_t k = 0; k < num_points; ++k)
    {
      result.push_back((-1 + 4 * k) * (-1 + 4 * k) * constant);
    }
  return result;
}
