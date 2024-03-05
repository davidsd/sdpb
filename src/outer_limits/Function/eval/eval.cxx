#include "chebyshev_clenshaw_recurrence.hxx"
#include "outer_limits/Function.hxx"

El::BigFloat
Function::eval(const El::BigFloat &epsilon, const El::BigFloat &infinity, const El::BigFloat &x) const
{
  // I use epsilon as a special value for when I want the limiting
  // value at zero.  If I use zero, then when I need to actually
  // evaluate it at zero, it gives the wrong answer.
  if(x==epsilon)
    {
      return epsilon_value;
    }
  else if(x == infinity)
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
