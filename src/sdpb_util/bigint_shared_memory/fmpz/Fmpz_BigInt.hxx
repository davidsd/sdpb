#pragma once

#include "sdpb_util/flint.hxx"

#include <El.hpp>

#include <boost/noncopyable.hpp>

// RAII wrapper for fmpz_t
// Only necessary members added
struct Fmpz_BigInt : private boost::noncopyable
{
  fmpz_t value;

  Fmpz_BigInt();
  ~Fmpz_BigInt();

  void from_BigFloat(const El::BigFloat &input);
  void to_BigFloat(El::BigFloat &output) const;
};
