#include "parse_Real.hxx"
#include "parse_many.hxx"

std::vector<Real> parse_vector(const boost::property_tree::ptree &tree)
{
  return parse_many("elt", parse_Real, tree);
}
