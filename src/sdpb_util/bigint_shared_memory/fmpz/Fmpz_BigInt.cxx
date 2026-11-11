#include "Fmpz_BigInt.hxx"

#include "fmpz_BigFloat_convert.hxx"

Fmpz_BigInt::Fmpz_BigInt()
{
  fmpz_init(value);
}
Fmpz_BigInt::~Fmpz_BigInt()
{
  fmpz_clear(value);
}

void Fmpz_BigInt::from_BigFloat(const El::BigFloat &input)
{
  BigFloat_to_fmpz_t(input, value);
}
void Fmpz_BigInt::to_BigFloat(El::BigFloat &output) const
{
  fmpz_t_to_BigFloat(value, output);
}
