#pragma once

#include "sdpb_util/assert.hxx"

#include <boost/optional.hpp>

#include <sstream>
#include <string>

// Converts string containing memory size in [kilo|mega|giga]bytes
// into number of bytes, e.g.:
// 100 or 100B -> 100
// 100K or 100KB -> 102400
// 100M or 100MB -> 104857600
// 100G or 100GB -> 107374182400
//
// Based on example from
// https://theboostcpplibraries.com/boost.propertytree#ex.propertytree_03
struct String_To_Bytes_Translator
{
  typedef std::string internal_type;
  typedef size_t external_type;

  static size_t from_string(const std::string &s)
  {
    std::istringstream iss(s);
    double number;
    std::string suffix;

    iss >> number >> suffix;

    size_t factor;
    if(suffix.empty() || suffix == "B")
      factor = 1;
    else if(suffix == "K" || suffix == "KB")
      factor = 1024;
    else if(suffix == "M" || suffix == "MB")
      factor = 1024 * 1024;
    else if(suffix == "G" || suffix == "GB")
      factor = 1024 * 1024 * 1024;
    else
      RUNTIME_ERROR("Cannot parse memory size: \"", s, "\"");

    return number * factor;
  }

  boost::optional<size_t> get_value(const std::string &s)
  {
    return from_string(s);
  }

  std::string put_value(size_t value) { return std::to_string(value); }
};
