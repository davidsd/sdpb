#pragma once

#include <El.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/noncopyable.hpp>
#include <boost/version.hpp>

#include <cmath>

using Boost_Float = boost::multiprecision::mpfr_float;

// Boost_Float's default precision for the calling thread, in decimal
// digits. Since Boost 1.77 each thread has its own default, copied from the
// global one (Boost_Float::default_precision) when the thread first uses
// Boost_Float; setting the global one later does not change it.
inline unsigned boost_float_thread_default_digits10()
{
#if BOOST_VERSION >= 107700
  return Boost_Float::thread_default_precision();
#else
  return Boost_Float::default_precision();
#endif
}
inline void set_boost_float_thread_default_digits10(const unsigned digits10)
{
#if BOOST_VERSION >= 107700
  Boost_Float::thread_default_precision(digits10);
#else
  Boost_Float::default_precision(digits10);
#endif
}

// Sets Boost_Float's default precision for the calling thread for a scope,
// in bits (as Environment::set_precision). Expressions with unrelated
// operands (strings, integers) are evaluated at the default precision, and
// on older Boost all of them are, so a computation that needs more bits
// than SDPB's precision sets it here.
struct Scoped_Boost_Float_Precision : boost::noncopyable
{
  explicit Scoped_Boost_Float_Precision(const mp_bitcnt_t precision)
      : old_digits10(boost_float_thread_default_digits10())
  {
    set_boost_float_thread_default_digits10(
      std::ceil(precision * std::log10(2.0)));
  }
  ~Scoped_Boost_Float_Precision()
  {
    set_boost_float_thread_default_digits10(old_digits10);
  }

private:
  unsigned old_digits10;
};

std::string to_string(const Boost_Float &boost_float);

Boost_Float to_Boost_Float(const El::BigFloat &alpha);

El::BigFloat to_BigFloat(const Boost_Float &value);

std::vector<El::BigFloat>
to_BigFloat_Vector(const std::vector<Boost_Float> &input);

std::vector<Boost_Float>
to_Boost_Float_Vector(const std::vector<El::BigFloat> &input);
