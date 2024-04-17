#pragma once

#include <istream>
#include <string>

#include <boost/algorithm/string.hpp>

enum Verbosity
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
    value = none;
  else if(token == "1" || token == "regular")
    value = regular;
  else if(token == "2" || token == "debug")
    value = debug;
  else if(token == "3" || token == "trace")
    value = trace;
  else
    in.setstate(std::ios_base::failbit);

  return in;
}
