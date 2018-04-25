#pragma once

#include <El.hpp>
#include <boost/property_tree/ptree.hpp>

std::vector<El::BigFloat>
parse_vector_elemental(const boost::property_tree::ptree &tree);
