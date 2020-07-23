#include "../Polynomial.hxx"

void load_poly(const int64_t &degree, const std::string &coefficient,
               Polynomial &poly)
{
  if(poly.coefficients.size() <= degree)
    {
      poly.coefficients.resize(degree + 1, El::BigFloat(0));
    }
  poly.coefficients[degree]=coefficient;
}
