#pragma once

#include "Damped_Rational_State.hxx"
#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <vector>
#include <string>
#include <stdexcept>

struct Positive_Matrix_With_Prefactor_State
{
  std::string name;
  bool inside = false, parsing_damped_rational = false,
    parsing_polynomials = false;
  Positive_Matrix_With_Prefactor value;

  Damped_Rational_State damped_rational_state;
  Vector_State<
    Vector_State<Vector_State<Vector_State<Number_State<El::BigFloat>>>>>
    polynomials_state;
  Positive_Matrix_With_Prefactor_State(const std::vector<std::string> &names,
                                       const size_t &offset)
      : name(names.at(offset)), damped_rational_state(names, offset + 1),
        polynomials_state(names, offset + 2)
  {}
  Positive_Matrix_With_Prefactor_State(
    const std::initializer_list<std::string> &names)
      : Positive_Matrix_With_Prefactor_State(names, 0)
  {}

  void json_key(const std::string &key);
  void json_string(const std::string &s);
  void json_start_array();
  void json_end_array();
  void json_start_object();
  void json_end_object();
};
