#include "parse_many.hxx"
#include "parse_Real.hxx"

std::vector<Real>
parse_vector(const boost::property_tree::ptree &tree)
{
  std::function<Real(const boost::property_tree::ptree &)> p (parse_Real);
  return parse_many("elt", p, tree);
}
