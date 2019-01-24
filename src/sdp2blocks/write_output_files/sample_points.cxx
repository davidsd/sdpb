#include <El.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfr.hpp>

#include <vector>

// Rescaled Laguerre
std::vector<El::BigFloat> sample_points(const size_t &num_points)
{
  std::vector<El::BigFloat> result;
  result.reserve(num_points);

  // Use boost::multiprecision because El::Bigfloat does not implement log().
  using boost_float = boost::multiprecision::mpfr_float;
  boost_float::default_precision(El::gmp::Precision());

  // Use raw mpfr to get π because boost::multiprecision crashes
  // before main() when trying to compute π for dynamic precision.
  // The precision has not been set, so episilon==0.  Computing π
  // requires 1/epsilon, which throws an exception.
  boost_float pi;
  mpfr_const_pi(pi.backend().data(),MPFR_RNDN);

  boost_float rho_crossing(3 - 2 * sqrt(boost_float(2))),
    constant_factor(-pi*pi / (64 * num_points * log(rho_crossing)));

  std::stringstream ss;
  ss.precision(El::gmp::Precision());
  ss << constant_factor;
  
  El::BigFloat constant;
  ss >> constant;

  std::cout << "constant: " << constant << "\n";
  
  for(size_t k = 0; k < num_points; ++k)
    {
      result.push_back((-1 + 4 * k) * (-1 + 4 * k) * constant);
    }
  return result;
}
