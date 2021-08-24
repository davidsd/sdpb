#pragma once

#include "Polynomial_Vector_Matrix_State.hxx"

using namespace std::string_literals;
class Input_Parser
{
public:
  std::string sdp_name = "sdp";
  bool inside_sdp = false, finished_objective = false,
       finished_polynomial_vector_matrices = false;

  Vector_State<Number_State<El::BigFloat>> objective_state;
  Vector_State<Polynomial_Vector_Matrix_State> polynomial_vector_matrices_state;

  Input_Parser(std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices,
               size_t &num_processed)
      : objective_state({"objective"s, "elt"s}),
        polynomial_vector_matrices_state(
          {"polynomialVectorMatrices"s, "polynomialVectorMatrix"s},
          polynomial_vector_matrices, num_processed)
  {}

  void on_start_element(const std::string &element_name);
  void on_end_element(const std::string &element_name);
  void on_characters(const xmlChar *characters, int length);
};
