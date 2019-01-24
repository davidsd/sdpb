#pragma once

#include <El.hpp>

#include <libxml2/libxml/parser.h>
#include <vector>
#include <sstream>
#include <stdexcept>

class Number_State
{
public:
  bool inside = false;
  // Need a string intermediate value because the parser may give the
  // element in chunks.  We need to concatenate them together
  // ourselves.
  std::stringstream string_value;
  El::BigFloat value;
  std::string name;

  Number_State(const std::string &Name) : name(Name) {}
  Number_State(const std::vector<std::string> &names, const size_t &offset)
      : Number_State(names.at(offset))
  {}
  Number_State() = delete;

  bool on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        throw std::runtime_error("Invalid input file.  Unexpected element '"
                                 + element_name + "' inside '" + name + "'");
      }
    else
      {
        inside = (element_name == name);
        if(inside)
          {
            string_value.str("");
            string_value.clear();
          }
      }
    return inside;
  }

  bool on_end_element(const std::string &)
  {
    bool result(inside);
    if(inside)
      {
        inside = false;
        string_value >> value;
      }
    return result;
  }

  bool on_characters(const xmlChar *characters, int length)
  {
    if(inside)
      {
        string_value.write(reinterpret_cast<const char *>(characters), length);
      }
    return inside;
  }
};
