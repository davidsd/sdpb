#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>

El::BigFloat power_prefactor(const Boost_Float &base, const El::BigFloat &x)
{
  // TODO: This is really, really inefficient
  // TODO: Need to add in constant term
  std::stringstream ss;
  set_stream_precision(ss);
  ss << x;
  Boost_Float x_mpfr(ss.str());
  return to_string(pow(base, x_mpfr));
}
