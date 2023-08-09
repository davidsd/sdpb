#pragma once

#include <El.hpp>
#include <boost/noncopyable.hpp>
#include <vector>

using Float = El::BigFloat;
using Float_Vector = std::vector<Float>;
using Float_Matrix = El::Matrix<Float>;

// Temporarily sets El::BigFloat precision
// and precision used in BigFloat comparison
struct Float_Binary_Precision : boost::noncopyable
{
  explicit Float_Binary_Precision(unsigned int binary_precision,
                                  unsigned int diff_precision);
  ~Float_Binary_Precision();

private:
  unsigned int old_binary_precision;
  unsigned int old_diff_precision;
};