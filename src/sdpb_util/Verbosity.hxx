#pragma once

#include "assert.hxx"

#include <istream>
#include <string>

#include <boost/algorithm/string.hpp>

enum class Verbosity
{
  none = 0,
  regular = 1,
  debug = 2,
  trace = 3,
};

inline std::istream &operator>>(std::istream &in, Verbosity &value)
{
  std::string token;
  in >> token;
  boost::algorithm::to_lower(token);

  if(token == "0" || token == "none")
    value = Verbosity::none;
  else if(token == "1" || token == "regular")
    value = Verbosity::regular;
  else if(token == "2" || token == "debug")
    value = Verbosity::debug;
  else if(token == "3" || token == "trace")
    value = Verbosity::trace;
  else
    in.setstate(std::ios_base::failbit);

  return in;
}

inline std::ostream &operator<<(std::ostream &os, const Verbosity verbosity)
{
  switch(verbosity)
    {
    case Verbosity::none: return os << "none";
    case Verbosity::regular: return os << "regular";
    case Verbosity::debug: return os << "debug";
    case Verbosity::trace: return os << "trace";
    default: THROW(std::out_of_range, "Block_File_Format=", (int)verbosity);
    }
  return os;
}
