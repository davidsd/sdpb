#pragma once

#include "Damped_Rational_State.hxx"
#include "../Positive_Matrix_With_Prefactor.hxx"

#include <libxml2/libxml/parser.h>
#include <vector>
#include <string>
#include <stdexcept>

struct Positive_Matrix_With_Prefactor_State
{
  std::string name;
  bool inside = false;
  Positive_Matrix_With_Prefactor value;

  Damped_Rational_State damped_rational_state;
  using Polynomial_State = Vector_State<Number_State>;
  // using Polynomial_State = Vector_State<Polynomial_Term_State>;
  Vector_State<Vector_State<Vector_State<Polynomial_State>>> polynomials_state;

  Positive_Matrix_With_Prefactor_State(const std::vector<std::string> &names,
                                       const size_t &offset)
      : name(names.at(offset)), damped_rational_state(names, offset + 1),
        polynomials_state(names, offset + 2)
  {}
  Positive_Matrix_With_Prefactor_State(
    const std::initializer_list<std::string> &names)
      : Positive_Matrix_With_Prefactor_State(names, 0)
  {}

  bool on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        if(element_name != "Symbol"
           && !damped_rational_state.on_start_element(element_name)
           && !polynomials_state.on_start_element(element_name))
          {
            throw std::runtime_error(
              "Invalid input file.  Expected 'Function' inside "
              "PositiveMatrixWithPrefactor, but found '"
              + element_name + "'");
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
        if(damped_rational_state.on_end_element(element_name))
          {
            using namespace std;
            swap(damped_rational_state.value, value.damped_rational);
          }
        else if(polynomials_state.on_end_element(element_name))
          {
            using namespace std;
            swap(value.polynomials,
                 polynomials_state.element_state.element_state.value);
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
        damped_rational_state.on_characters(characters, length)
          || polynomials_state.on_characters(characters, length);
      }
    return inside;
  }
};
