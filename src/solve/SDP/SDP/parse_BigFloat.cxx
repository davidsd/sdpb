#include <El.hpp>

#include <boost/property_tree/ptree.hpp>

El::BigFloat parse_BigFloat(const boost::property_tree::ptree &tree)
{
  return El::BigFloat(tree.data(), 10);
}
