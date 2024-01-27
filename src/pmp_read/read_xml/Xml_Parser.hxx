#pragma once

#include "Xml_Polynomial_Vector_Matrix_State.hxx"

using namespace std::string_literals;
class Xml_Parser
{
public:
  std::string sdp_name = "sdp";
  bool inside_sdp = false, finished_objective = false,
       finished_polynomial_vector_matrices = false;

  Vector_State<Number_State<El::BigFloat>> objective_state;
  Vector_State<Xml_Polynomial_Vector_Matrix_State> polynomial_vector_matrices_state;

  Xml_Parser(const std::function<void(Xml_Polynomial_Vector_Matrix_State& matrix_state)> &process_matrix)
      : objective_state({"objective"s, "elt"s}),
        polynomial_vector_matrices_state(
          {"polynomialVectorMatrices"s, "polynomialVectorMatrix"s},
          process_matrix)
  {}

  void on_start_element(const std::string &element_name);
  void on_end_element(const std::string &element_name);
  void on_characters(const xmlChar *characters, int length);
};
