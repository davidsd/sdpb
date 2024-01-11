#pragma once

#include "Damped_Rational_State.hxx"
#include "pmp/Polynomial_Vector_Matrix.hxx"

#include <memory>
#include <vector>
#include <string>

struct Json_Polynomial_Vector_Matrix_State
{
  std::string name;
  bool inside = false, parsing_damped_rational = false,
       parsing_polynomials = false;
  std::unique_ptr<Polynomial_Vector_Matrix> value;

  Damped_Rational_State damped_rational_state;
  Vector_State<
    Vector_State<Vector_State<Vector_State<Number_State<El::BigFloat>>>>>
    polynomials_state;
  Json_Polynomial_Vector_Matrix_State(const std::vector<std::string> &names,
                                      const size_t &offset)
      : name(names.at(offset)),
        damped_rational_state(names, offset + 1),
        polynomials_state(names, offset + 2)
  {}
  Json_Polynomial_Vector_Matrix_State(
    const std::initializer_list<std::string> &names)
      : Json_Polynomial_Vector_Matrix_State(names, 0)
  {}

  void json_key(const std::string &key);
  void json_string(const std::string &s);
  void json_start_array();
  void json_end_array();
  void json_start_object();
  void json_end_object();
};
