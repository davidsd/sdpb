#include "../Polynomial_Vector_Matrix.hxx"
#include "parse_Real.hxx"
#include "parse_many.hxx"
#include "parse_vector.hxx"

namespace
{
  Polynomial parse_polynomial(const boost::property_tree::ptree &tree)
  {
    Polynomial result;
    result.coefficients = parse_many("coeff", parse_Real, tree);
    return result;
  }

  std::vector<Polynomial>
  parse_polynomial_vector(const boost::property_tree::ptree &tree)
  {
    return parse_many("polynomial", parse_polynomial, tree);
  }
}

Polynomial_Vector_Matrix
parse_polynomial_vector_matrix(const boost::property_tree::ptree &tree)
{
  Polynomial_Vector_Matrix result;
  result.rows = tree.get<int64_t>("rows");
  result.cols = tree.get<int64_t>("cols");

  result.elements = parse_many("polynomialVector", parse_polynomial_vector,
                               tree.get_child("elements"));
  result.samplePoints = parse_vector(tree.get_child("samplePoints"));
  result.sampleScalings = parse_vector(tree.get_child("sampleScalings"));
  result.bilinearBasis = parse_many("polynomial", parse_polynomial,
                                    tree.get_child("bilinearBasis"));
  return result;
}
