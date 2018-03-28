#pragma once

#include "../types.hxx"

#include <boost/property_tree/ptree.hpp>

#include <vector>

std::vector<Real>
parse_vector(const boost::property_tree::ptree &tree);
