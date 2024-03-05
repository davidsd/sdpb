#pragma once

#include <istream>
#include <string>

enum Verbosity
{
  none = 0,
  regular = 1,
  debug = 2
};

inline std::istream &operator>>(std::istream &in, Verbosity &value)
{
  std::string token;
  in >> token;

  if(token == "0" || token == "none")
    value = none;
  else if(token == "1" || token == "regular")
    value = regular;
  else if(token == "2" || token == "debug")
    value = debug;
  else
    in.setstate(std::ios_base::failbit);

  return in;
}
