#include "../../Boost_Float.hxx"
#include "../../../Polynomial.hxx"

std::vector<Polynomial>
bilinear_form(const Damped_Rational &damped_rational, const int64_t &m)
{
  std::vector<Boost_Float> sorted_poles(damped_rational.poles);
  std::sort(sorted_poles.begin(), sorted_poles.end());
  for(auto pole(sorted_poles.begin()); pole != sorted_poles.end();)
    {
      El::BigFloat &p(*pole);
      auto range(std::equal_range(pole,sorted_poles.end(),p));
      auto l(std::distance(range.first,range.second));

      
      do
        {
          ++pole;
        }
      while(Abs(p-*pole)<1.0e-2);
    }
}
