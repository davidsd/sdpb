#pragma once

#include "sdpb_util/assert.hxx"

#include <iostream>
#include <string>

enum Block_File_Format
{
  bin,
  json
};

inline std::istream &operator>>(std::istream &in, Block_File_Format &format)
{
  std::string token;
  in >> token;
  if(token == "bin")
    format = bin;
  else if(token == "json")
    format = json;
  else
    in.setstate(std::ios_base::failbit);
  return in;
}

inline std::ostream &
operator<<(std::ostream &out, const Block_File_Format &format)
{
  if(format == bin)
    out << "bin";
  else if(format == json)
    out << "json";
  else
    THROW(std::out_of_range, "Block_File_Format=", std::to_string(format));
  return out;
}
