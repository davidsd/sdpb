#pragma once

#include "Positive_Matrix_With_Prefactor_State.hxx"

using namespace std::string_literals;
class Input_Parser
{
public:
  std::string expression_name = "Expression", sdp_name = "Function";
  bool inside_expression = false, inside_sdp = false,
       finished_objective = false, finished_normalization = false;

  Vector_State<Number_State> objective_state, normalization_state;
  // Vector_State<Positive_Matrix_With_Prefactor_State>
  //   positive_matrix_with_prefactor_state;

  Input_Parser()
      : objective_state({"Function"s, "Number"s}),
        normalization_state({"Function"s, "Number"s})
      // ,
      //   positive_matrix_with_prefactor_state(
      //     {"Function"s, "Function"s, "Function"s, "Function"s, "Function"s,
      //      "Function"s, "Function"s, "Function"s})
  {}

  void on_start_element(const std::string &element_name);
  void on_end_element(const std::string &element_name);
  void on_characters(const xmlChar *characters, int length);
};
