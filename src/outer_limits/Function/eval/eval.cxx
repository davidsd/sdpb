#include "chebyshev_clenshaw_recurrence.hxx"
#include "../../Function.hxx"

El::BigFloat
Function::eval(const El::BigFloat &infinity, const El::BigFloat &x) const
{
  if(x == infinity)
    {
      // std::cout << "eval inf: "
      //           << infinity_value << "\n";
      return infinity_value;
    }
  else
    {
      // std::cout << "eval: "
      //           << max_delta << " "
      //           << chebyshev_values.size() << " "
      //           << x << " "
      //           << chebyshev_clenshaw_recurrence(chebyshev_values.data(),
      //                                            chebyshev_values.size(),
      //                                            El::BigFloat(0.0), max_delta, x) << "\n";
      return chebyshev_clenshaw_recurrence(chebyshev_coeffs.data(),
                                           chebyshev_coeffs.size(),
                                           El::BigFloat(0.0), max_delta, x);
    }
}
