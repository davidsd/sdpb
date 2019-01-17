#pragma once

#include "../../Vector_State.hxx"
#include "../../Number_State.hxx"
#include "../Damped_Rational.hxx"

#include <libxml2/libxml/parser.h>
#include <vector>
#include <string>
#include <stdexcept>

using namespace std::string_literals;
struct Damped_Rational_State
{
  std::string name;
  bool inside = false, finished_constant = false, finished_base = false;
  Damped_Rational value;
  Number_State constant_state, base_state;
  Vector_State<Number_State> polynomial_state;

  Damped_Rational_State(const std::vector<std::string> &names,
                        const size_t &offset)
      : name(names.at(offset)), constant_state("Number"s),
        base_state("Number"s), polynomial_state({"Function"s, "Number"s})
  {}
  Damped_Rational_State(const std::initializer_list<std::string> &names)
      : Damped_Rational_State(names, 0)
  {}

  bool on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        if(element_name != "Symbol"
           && (finished_constant
               || !constant_state.on_start_element(element_name))
           && !polynomial_state.on_start_element(element_name)
           && (finished_base || !base_state.on_start_element(element_name)))
          {
            throw std::runtime_error(
              "Invalid input file.  Expected '" + constant_state.name
              + "' or '" + polynomial_state.name
              + "' inside DampedRational, but found '" + element_name + "'");
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
        if(element_name != "Symbol")
          {
            if(constant_state.on_end_element(element_name))
              {
                finished_constant = !constant_state.inside;
                if(finished_constant)
                  {
                    using namespace std;
                    swap(value.constant, constant_state.value);
                  }
              }
            else if(polynomial_state.on_end_element(element_name))
              {
                if(!polynomial_state.inside)
                  {
                    using namespace std;
                    swap(value.poles, polynomial_state.value);
                  }
              }
            else if(base_state.on_end_element(element_name))
              {
                finished_base = !base_state.inside;
                if(finished_base)
                  {
                    using namespace std;
                    swap(value.base, base_state.value);
                  }
              }
            else
              {
                inside = (element_name != name);
              }
          }
      }
    return result;
  }

  bool on_characters(const xmlChar *characters, int length)
  {
    if(inside)
      {
        constant_state.on_characters(characters, length)
          || polynomial_state.on_characters(characters, length)
          || base_state.on_characters(characters, length);
      }
    return inside;
  }
};
