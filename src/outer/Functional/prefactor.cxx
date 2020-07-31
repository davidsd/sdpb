#include "../Functional.hxx"
#include "../../set_stream_precision.hxx"

#include "../../sdp2input/Boost_Float.hxx"

El::BigFloat Functional::prefactor(const El::BigFloat &x) const
{
  std::stringstream ss;
  set_stream_precision(ss);
  ss << x;
  Boost_Float x_mpfr(ss.str());
  return to_string(pow(prefactor_power,x_mpfr));
}
