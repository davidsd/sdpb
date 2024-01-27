#pragma once

#include "assert.hxx"

#include <El.hpp>

#include <libxml2/libxml/parser.h>
#include <vector>
#include <sstream>
#include <stdexcept>

template <typename Float_Type> class Number_State
{
public:
  bool inside = false;
  // Need a string intermediate value because the parser may give the
  // element in chunks.  We need to concatenate them together
  // ourselves.
  std::stringstream string_value;
  Float_Type value;
  std::string name;

  Number_State(const std::string &Name) : name(Name) {}
  Number_State(const std::vector<std::string> &names, const size_t &offset)
      : Number_State(names.at(offset))
  {}
  Number_State() = delete;

  // XML Functions
  bool xml_on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        RUNTIME_ERROR("Invalid input file. Unexpected element '", element_name,
                      "' inside '", name, "'");
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

  bool xml_on_end_element(const std::string &)
  {
    bool result(inside);
    if(inside)
      {
        inside = false;
        try
          {
            // GMP does not have inf or nan, so we approximate inf
            // with max double.
            using namespace std::string_literals;
            if(string_value.str() == "inf"s)
              {
                value = Float_Type(std::numeric_limits<double>::max());
              }
            else
              {
                value = Float_Type(string_value.str());
              }
          }
        catch(...)
          {
            RUNTIME_ERROR("Invalid NNNN number: '", string_value.str(), "'");
          }
      }
    return result;
  }

  bool xml_on_characters(const xmlChar *characters, int length)
  {
    if(inside)
      {
        string_value.write(reinterpret_cast<const char *>(characters), length);
      }
    return inside;
  }

  // JSON Functions
  void json_key(const std::string &key)
  {
    RUNTIME_ERROR("Invalid input file. Found the key '", key,
                  "' when expecting a number.");
  }

  void json_string(const std::string &s)
  {
    try
      {
        // GMP does not have inf or nan, so we approximate inf
        // with max double.
        using namespace std::string_literals;
        if(s == "inf"s)
          {
            value = Float_Type(std::numeric_limits<double>::max());
          }
        else
          {
            value = Float_Type(s);
          }
      }
    catch(...)
      {
        RUNTIME_ERROR("Invalid number: '", s, "'");
      }
  }

  void json_start_array()
  {
    RUNTIME_ERROR(
      "Invalid input file. Found an array when expecting a number.");
  }

  void json_end_array()
  {
    RUNTIME_ERROR(
      "Invalid input file. Found an array end when parsing a number.");
  }

  void json_start_object()
  {
    RUNTIME_ERROR(
      "Invalid input file. Found an object when expecting a number.");
  }

  void json_end_object()
  {
    RUNTIME_ERROR(
      "Invalid input file. Found an object end when parsing a number.");
  }
};
