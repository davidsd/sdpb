#pragma once

#include <El.hpp>

inline void fmpz_t_to_BigFloat(const fmpz_t input, El::BigFloat &output)
{
  fmpz_get_mpf(output.gmp_float.get_mpf_t(), input);
}
inline void BigFloat_to_fmpz_t(const El::BigFloat &input, fmpz_t output)
{
  fmpz_set_mpf(output, input.gmp_float.get_mpf_t());
}
