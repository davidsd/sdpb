#include "../Functional.hxx"
#include "../../set_stream_precision.hxx"

#include "../../sdp2input/Boost_Float.hxx"

El::BigFloat Functional::prefactor(const El::BigFloat &x)
{
  if(!has_prefactor)
    {
      return El::BigFloat(1);
    }

  std::stringstream ss;
  set_stream_precision(ss);
  ss << x;
  Boost_Float x_mpfr(ss.str());
  return to_string(pow(4 * (3 - 2 * sqrt(Boost_Float(2.0))), x_mpfr));
}
