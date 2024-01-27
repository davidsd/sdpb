#pragma once

#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"

#include <El.hpp>
#include <vector>

struct Damped_Rational
{
  Boost_Float constant, base;
  std::vector<Boost_Float> poles;

  bool is_constant() const
  {
    return poles.empty() && base==1;
  }
};

inline void swap(Damped_Rational &a, Damped_Rational &b)
{
  using namespace std;
  swap(a.constant, b.constant);
  swap(a.base, b.base);
  swap(a.poles, b.poles);
}

inline std::ostream & operator<<(std::ostream &os, const Damped_Rational &damped)
{
  os << "{\n  \"constant\": \""
     << damped.constant << "\",\n"
     << "  \"base\": \""
     << damped.base << "\",\n"
     << "  \"poles\": "
     << damped.poles << "\n}";
  return os;
}
