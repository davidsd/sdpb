#pragma once

#include <libxml2/libxml/parser.h>
#include <vector>
#include <string>
#include <stdexcept>

template <typename T> class Vector_State
{
public:
  std::string name;
  bool inside = false;
  T element_state;
  std::vector<decltype(element_state.value)> value;

  Vector_State(const std::vector<std::string> &names, const size_t &offset)
      : name(names.at(offset)), element_state(names, offset + 1)
  {}
  template <typename U>
  Vector_State(const std::vector<std::string> &names, const size_t &offset,
               U &u)
      : name(names.at(offset)), element_state(names, offset + 1, u)
  {}
  template <typename U, typename V>
  Vector_State(const std::vector<std::string> &names, const size_t &offset,
               U &u, V &v)
      : name(names.at(offset)), element_state(names, offset + 1, u, v)
  {}

  Vector_State(const std::initializer_list<std::string> &names)
      : Vector_State(names, 0)
  {}
  template <typename U>
  Vector_State(const std::initializer_list<std::string> &names, U &u)
      : Vector_State(names, 0, u)
  {}
  template <typename U, typename V>
  Vector_State(const std::initializer_list<std::string> &names, U &u, V &v)
      : Vector_State(names, 0, u, v)
  {}

  // XML Functions
  bool xml_on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        if(!element_state.xml_on_start_element(element_name))
          {
            throw std::runtime_error("Invalid input file.  Expected '"
                                     + element_state.name + "' inside Vector '"
                                     + name + "', but found '" + element_name
                                     + "'");
          }
      }
    else
      {
        inside = (element_name == name);
        if(inside)
          {
            value.clear();
          }
      }
    return inside;
  }

  bool xml_on_end_element(const std::string &element_name)
  {
    bool result(inside);
    if(inside)
      {
        if(element_state.xml_on_end_element(element_name))
          {
            if(!element_state.inside)
              {
                value.emplace_back();
                std::swap(value.back(), element_state.value);
              }
          }
        else
          {
            inside = (element_name != name);
          }
      }
    return result;
  }

  bool xml_on_characters(const xmlChar *characters, int length)
  {
    if(inside)
      {
        element_state.xml_on_characters(characters, length);
      }
    return inside;
  }

  // JSON Functions
  void json_key(const std::string &key) { element_state.json_key(key); }

  void json_string(const std::string &s)
  {
    element_state.json_string(s);
    if(!element_state.inside)
      {
        value.emplace_back();
        std::swap(value.back(), element_state.value);
      }
  }

  void json_start_array()
  {
    if(inside)
      {
        element_state.json_start_array();
      }
    else
      {
        inside = true;
        value.clear();
      }
  }

  void json_end_array()
  {
    if(element_state.inside)
      {
        element_state.json_end_array();
        if(!element_state.inside)
          {
            value.emplace_back();
            std::swap(value.back(), element_state.value);
          }
      }
    else
      {
        inside = false;
      }
  }

  void json_start_object() { element_state.json_start_object(); }

  void json_end_object()
  {
    element_state.json_end_object();
    if(!element_state.inside)
      {
        value.emplace_back();
        std::swap(value.back(), element_state.value);
      }
  }
};
