#include "parse_BigFloat.hxx"
#include "parse_many.hxx"

std::vector<El::BigFloat> parse_vector(const boost::property_tree::ptree &tree)
{
  return parse_many("elt", parse_BigFloat, tree);
}
