#pragma once

#include "Boost_Float.hxx"
#include <El.hpp>
#include <vector>

struct Damped_Rational
{
  Boost_Float constant, base;
  std::vector<Boost_Float> poles;
  std::string variable;
};

inline void swap(Damped_Rational &a, Damped_Rational &b)
{
  using namespace std;
  swap(a.constant, b.constant);
  swap(a.base, b.base);
  swap(a.poles, b.poles);
  swap(a.variable, b.variable);
}
