#include "../../PolynomialVectorMatrix.hxx"
#include "../parse_Real.hxx"
#include "../parse_many.hxx"
#include "../parse_vector.hxx"

namespace
{
  Polynomial parse_polynomial(const boost::property_tree::ptree &tree)
  {
    Polynomial result;
    std::function<Real(const boost::property_tree::ptree &)> p(parse_Real);
    result.coefficients = parse_many("coeff", p, tree);
    return result;
  }

  std::vector<Polynomial>
  parse_polynomial_vector(const boost::property_tree::ptree &tree)
  {
    std::function<Polynomial(const boost::property_tree::ptree &)> p(
    parse_polynomial);
    return parse_many("polynomial", p, tree);
  }
}

PolynomialVectorMatrix
parse_polynomial_vector_matrix(const boost::property_tree::ptree &tree)
{
  PolynomialVectorMatrix result;
  result.rows = tree.get<int64_t>("rows");
  result.cols = tree.get<int64_t>("cols");
  std::function<std::vector<Polynomial>(const boost::property_tree::ptree &)>
  pv(parse_polynomial_vector);

  result.elements
  = parse_many("polynomialVector", pv, tree.get_child("elements"));
  result.samplePoints = parse_vector(tree.get_child("samplePoints"));
  result.sampleScalings = parse_vector(tree.get_child("sampleScalings"));
  std::function<Polynomial(const boost::property_tree::ptree &)> pp(
  parse_polynomial);
  result.bilinearBasis
  = parse_many("polynomial", pp, tree.get_child("bilinearBasis"));
  return result;
}
