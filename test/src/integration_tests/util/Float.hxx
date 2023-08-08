#pragma once

#include <boost/multiprecision/mpfr.hpp>
#include <vector>

using Float = boost::multiprecision::mpfr_float;
using Float_Vector = std::vector<Float>;
using Float_Matrix = std::vector<Float_Vector>;

// Temporarily sets MPFR precision
struct Float_Binary_Precision : boost::noncopyable
{
  explicit Float_Binary_Precision(unsigned int binary_precision);
  ~Float_Binary_Precision();

private:
  unsigned int old_decimal_precision;
};