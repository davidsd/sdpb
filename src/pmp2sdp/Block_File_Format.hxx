#pragma once

#include "sdpb_util/assert.hxx"

#include <iostream>
#include <string>

enum class Block_File_Format
{
  bin,
  json
};

inline std::istream &operator>>(std::istream &in, Block_File_Format &format)
{
  std::string token;
  in >> token;
  if(token == "bin")
    format = Block_File_Format::bin;
  else if(token == "json")
    format = Block_File_Format::json;
  else
    in.setstate(std::ios_base::failbit);
  return in;
}

inline std::ostream &
operator<<(std::ostream &out, const Block_File_Format &format)
{
  if(format == Block_File_Format::bin)
    out << "bin";
  else if(format == Block_File_Format::json)
    out << "json";
  else
    THROW(std::out_of_range, "Block_File_Format=", (int)format);
  return out;
}
