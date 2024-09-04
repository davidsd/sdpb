#include "Boost_Float.hxx"

#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <boost/multiprecision/mpfr.hpp>

#include <mpf2mpfr.h>

#include <sstream>

std::string to_string(const Boost_Float &boost_float)
{
  std::stringstream ss;
  set_stream_precision(ss);
  ss << boost_float;
  return ss.str();
}
Boost_Float to_Boost_Float(const El::BigFloat &alpha)
{
  Boost_Float result;
  mpfr_t &mpfr_alpha = result.backend().data();
  mpfr_set_prec(mpfr_alpha, alpha.gmp_float.get_prec());
  mpfr_set_f(mpfr_alpha, alpha.gmp_float.get_mpf_t(), MPFR_RNDN);
  return result;
}
El::BigFloat to_BigFloat(const Boost_Float &value)
{
  El::BigFloat result;
  mpfr_get_f(result.gmp_float.get_mpf_t(), value.backend().data(), MPFR_RNDN);
  return result;
}
std::vector<El::BigFloat>
to_BigFloat_Vector(const std::vector<Boost_Float> &input)
{
  std::vector<El::BigFloat> output;
  output.reserve(input.size());
  for(const auto &x : input)
    {
      output.push_back(to_BigFloat(x));
    }

  return output;
}
std::vector<Boost_Float>
to_Boost_Float_Vector(const std::vector<El::BigFloat> &input)
{
  std::vector<Boost_Float> output;
  output.reserve(input.size());
  for(const auto &x : input)
    {
      output.push_back(to_Boost_Float(x));
    }

  return output;
}