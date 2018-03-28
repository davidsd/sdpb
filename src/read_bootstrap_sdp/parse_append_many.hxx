#pragma once

#include "../types.hxx"

#include <boost/property_tree/ptree.hpp>

#include <vector>

template <class T>
std::vector<T> parse_append_many(const std::string &name,
                                 T (*parse_function)(const boost::property_tree::ptree &),
                                 const boost::property_tree::ptree &tree,
                                 std::vector<T> &result)
{
  for(auto &child : tree)
    {
      if(child.first == name)
        { result.emplace_back(parse_function(child.second)); } }
  return result;
}
