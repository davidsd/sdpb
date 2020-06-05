#pragma once

#include "Damped_Rational_State.hxx"
#include "Polynomial_State.hxx"
#include "../Positive_Matrix_With_Prefactor.hxx"

#include <libxml2/libxml/parser.h>
#include <vector>
#include <string>
#include <stdexcept>

struct Positive_Matrix_With_Prefactor_State
{
  std::string name;
  bool inside = false;
  // TODO: Remove duplication
  // For XML
  bool finished_damped_rational = false;
  // For JSON
  bool parsing_damped_rational = false, parsing_polynomials = false;
  Positive_Matrix_With_Prefactor value;

  Damped_Rational_State damped_rational_state;
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

  bool xml_on_start_element(const std::string &element_name);
  bool xml_on_end_element(const std::string &element_name);
  bool xml_on_characters(const xmlChar *characters, int length);

  void json_key(const std::string &key);
  void json_string(const std::string &s);
  void json_start_array();
  void json_end_array();
  void json_start_object();
  void json_end_object();
};
