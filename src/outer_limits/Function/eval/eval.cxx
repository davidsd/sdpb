#include "chebyshev_clenshaw_recurrence.hxx"
#include "../../Function.hxx"

El::BigFloat
Function::eval(const El::BigFloat &infinity, const El::BigFloat &x) const
{
  if(x == infinity)
    {
      return infinity_value;
    }
  else
    {
      return chebyshev_clenshaw_recurrence(chebyshev_coeffs.data(),
                                           chebyshev_coeffs.size(),
                                           El::BigFloat(0.0), max_delta, x);
    }
}
