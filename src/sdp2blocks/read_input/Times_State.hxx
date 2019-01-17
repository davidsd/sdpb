#pragma once

#include "Power_State.hxx"

struct Times_State
{
  std::string name;
  bool inside = false;
  std::pair<size_t, El::BigFloat> value;

  Number_State number_state;
  Power_State power_state;

  Times_State(const std::string &Name)
      : name(Name), number_state("Number"), power_state("Function")
  {}

  bool on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        // Have to check power_state before number_state, because inside
        // power_state there is a valid Number.
        if(element_name != "Symbol"
           && (power_state.inside
               || !number_state.on_start_element(element_name))
           // && (!power_state.inside
           //     && !number_state.on_start_element(element_name))
           && !power_state.on_start_element(element_name))
          {
            throw std::runtime_error(
              "Invalid input file.  Unexpected element '" + element_name
              + "' inside a Times element.");
          }
      }
    else
      {
        inside = (element_name == name);
        if(inside)
          {
            value.first = 1;
          }
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
                using namespace std;
                swap(value.second, number_state.value);
              }
          }
        else if(power_state.on_end_element(element_name))
          {
            if(!power_state.inside)
              {
                value.first = power_state.value;
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
        number_state.on_characters(characters, length)
          || power_state.on_characters(characters, length);
      }
    return inside;
  }
};
