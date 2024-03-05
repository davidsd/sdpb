#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>

El::BigFloat power_prefactor(const Boost_Float &base, const El::BigFloat &x)
{
  // TODO: Need to add in constant term
  return to_BigFloat(pow(base, to_Boost_Float(x)));
}
