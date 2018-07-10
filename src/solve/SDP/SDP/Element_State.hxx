#pragma once

#include <libxml++/libxml++.h>

#include <string>

class Element_State
{
public:
  bool inside = false;
  // Need a string intermediate value because the parser may give the
  // element in chunks.  We need to concatenate them together
  // ourselves.
  std::string string_value;
  El::BigFloat value;
  std::string name;

  Element_State(const std::string &Name) : name(Name) {}
  Element_State(const std::vector<std::string> &names, const size_t &offset)
      : Element_State(names.at(offset))
  {}
  Element_State() = delete;

  bool on_start_element(const Glib::ustring &element_name,
                        const xmlpp::SaxParser::AttributeList &)
  {
    if(inside)
      {
        throw std::runtime_error("Invalid input file.  Unexpected element '"
                                 + element_name + "' inside '" + name + "'");
      }
    else if(element_name == name)
      {
        inside = true;
        string_value.clear();
      }
    return inside;
  }

  bool on_end_element(const Glib::ustring &element_name)
  {
    bool result(false);
    if(inside)
      {
        if(element_name == name)
          {
            inside = false;
            result = true;
            value = El::BigFloat(string_value, 10);
          }
      }
    return result;
  }

  bool on_characters(const Glib::ustring &characters)
  {
    if(inside)
      {
        string_value += characters.data();
      }
    return inside;
  }
};
