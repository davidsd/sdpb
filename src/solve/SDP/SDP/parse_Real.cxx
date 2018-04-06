#include "../types.hxx"

#include <boost/property_tree/ptree.hpp>

Real parse_Real(const boost::property_tree::ptree &tree)
{
  return Real(tree.data());
}
