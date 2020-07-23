#pragma once

#include "../Polynomial.hxx"
#include "../ostream_vector.hxx"

struct Block
{
  std::vector<Polynomial> polys;
  std::vector<El::BigFloat> poles;

  void assign_poly(const size_t &poly_index, const size_t &degree_index,
                   const std::string &coefficient_string)

  {
    if(polys.size() <= poly_index)
      {
        polys.resize(poly_index + 1);
      }
    if(polys[poly_index].coefficients.size() <= degree_index)
      {
        polys[poly_index].coefficients.resize(degree_index + 1,
                                              El::BigFloat(0));
      }
    polys[poly_index].coefficients[degree_index] = coefficient_string;
  }

  El::BigFloat pole_prefactor(const El::BigFloat &x)
  {
    El::BigFloat pole_product(1);
    for(auto &pole : poles)
      {
        pole_product *= pole + x;
      }
    return 1 / pole_product;
  }

  El::BigFloat eval_weighted(const El::BigFloat &x,
                             const std::vector<El::BigFloat> &weights)
  {
    if(weights.size() != polys.size())
      {
        throw std::runtime_error("INTERNAL ERROR mismatch");
      }

    El::BigFloat result(0);
    for(size_t index(0); index != weights.size(); ++index)
      {
        result += weights[index] * polys[index](x);
      }
    if(!poles.empty())
      {
        result *= pole_prefactor(x);
      }
    return result;
  }
};

inline std::ostream & operator<<(std::ostream &os, const Block &block)
{
  os << block.polys;
  return os;
}
