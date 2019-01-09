#pragma once

#include "../../Vector_State.hxx"
#include "../../Element_State.hxx"

using namespace std::string_literals;
class Input_Parser
{
public:
  std::string expression_name = "Expression",
    sdp_name = "Function";
  bool inside_expression = false, inside_sdp = false,
       finished_objective = false, finished_normalization = false;

  Vector_State<Element_State> objective_state, normalization_state;

  Input_Parser()
      : objective_state({"Function"s, "Number"s}),
        normalization_state({"Function"s, "Number"s})
  {}

  void on_start_element(const std::string &element_name);
  void on_end_element(const std::string &element_name);
  void on_characters(const xmlChar *characters, int length);
};
