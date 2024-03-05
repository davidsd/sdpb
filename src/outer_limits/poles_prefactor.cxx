#include "sdpb_util/Boost_Float.hxx"

#include <El.hpp>

#include <vector>

El::BigFloat
poles_prefactor(const std::vector<Boost_Float> &poles, const El::BigFloat &x)
{
  El::BigFloat pole_product(1);
  for(auto &pole : poles)
    {
      pole_product *= x - to_BigFloat(pole);
    }
  return 1 / pole_product;
}
