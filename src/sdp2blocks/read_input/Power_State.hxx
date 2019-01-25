#pragma once

#include "../../Number_State.hxx"

struct Power_State
{
  std::string name;
  bool inside = false;
  size_t value;

  Number_State<El::BigFloat> number_state;

  Power_State(const std::string &Name) : name(Name), number_state("Number"s) {}
  Power_State() = delete;

  bool on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        if(element_name != "Symbol"
           && !number_state.on_start_element(element_name))
          {
            throw std::runtime_error(
              "Invalid input file.  Unexpected element '" + element_name
              + "' inside a Power element.");
          }
      }
    else
      {
        inside = (element_name == name);
      }
    return inside;
  }

  bool on_end_element(const std::string &element_name)
  {
    bool result(inside);
    if(inside)
      {
        if(number_state.on_end_element(element_name))
          {
            if(!number_state.inside)
              {
                number_state.string_value.seekg(0);
                number_state.string_value >> value;
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
        number_state.on_characters(characters, length);
      }
    return inside;
  }
};
