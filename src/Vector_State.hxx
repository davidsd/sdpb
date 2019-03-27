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
  template <typename U, typename V, typename W>
  Vector_State(const std::vector<std::string> &names, const size_t &offset,
               U &u, V &v, W &w)
      : name(names.at(offset)), element_state(names, offset + 1, u, v, w)
  {}

  Vector_State(const std::initializer_list<std::string> &names)
      : Vector_State(names, 0)
  {}
  template <typename U, typename V, typename W>
  Vector_State(const std::initializer_list<std::string> &names, U &u, V &v,
               W &w)
      : Vector_State(names, 0, u, v, w)
  {}

  bool on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        if(!element_state.on_start_element(element_name))
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

  bool on_end_element(const std::string &element_name)
  {
    bool result(inside);
    if(inside)
      {
        if(element_state.on_end_element(element_name))
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

  bool on_characters(const xmlChar *characters, int length)
  {
    if(inside)
      {
        element_state.on_characters(characters, length);
      }
    return inside;
  }
};
