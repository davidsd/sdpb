#pragma once

#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"

#include <El.hpp>
#include <vector>

struct Damped_Rational
{
  Boost_Float constant, base;
  std::vector<Boost_Float> poles;

  bool is_constant() const { return poles.empty() && base == 1; }

  // Evaluate at x, optionally regularize if x is close to the pole
  [[nodiscard]] Boost_Float
  evaluate(const Boost_Float &x,
           const Boost_Float &min_pole_distance = 0) const
  {
    const Boost_Float numerator = constant * pow(base, x);
    Boost_Float denominator = 1;
    for(auto &pole : poles)
      {
        Boost_Float delta = x - pole;

        // Regularize value at small denominator
        if(abs(delta) < min_pole_distance)
          {
            delta = min_pole_distance;
            if(sign(delta) < 0)
              delta = -delta;
          }

        denominator *= delta;
      }
    return numerator / denominator;
  }
};

inline void swap(Damped_Rational &a, Damped_Rational &b)
{
  using namespace std;
  swap(a.constant, b.constant);
  swap(a.base, b.base);
  swap(a.poles, b.poles);
}

inline std::ostream &
operator<<(std::ostream &os, const Damped_Rational &damped)
{
  os << "{\n  \"constant\": \"" << damped.constant << "\",\n"
     << "  \"base\": \"" << damped.base << "\",\n"
     << "  \"poles\": " << damped.poles << "\n}";
  return os;
}
