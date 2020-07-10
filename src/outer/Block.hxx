#pragma once

#include "../Polynomial.hxx"

struct Block
{
  std::vector<Polynomial> p;
  std::vector<El::BigFloat> poles;

  El::BigFloat
  eval(const El::BigFloat &x, const std::vector<El::BigFloat> &optimal)
  {
    if(optimal.size() != p.size())
      {
        throw std::runtime_error("INTERNAL ERROR mismatch");
      }
    El::BigFloat pole_product(1);
    for(auto &pole : poles)
      {
        pole_product *= pole + x;
      }
    const El::BigFloat pole_inverse(1 / pole_product);

    El::BigFloat result(0);
    for(size_t index(0); index != optimal.size(); ++index)
      {
        result+= optimal[index] * p[index](x) * pole_inverse;
      }
    return result;
  }
};
