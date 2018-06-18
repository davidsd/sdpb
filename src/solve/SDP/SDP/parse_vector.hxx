#pragma once

#include <El.hpp>
#include <boost/property_tree/ptree.hpp>

std::vector<El::BigFloat>
parse_vector(const boost::property_tree::ptree &tree);
