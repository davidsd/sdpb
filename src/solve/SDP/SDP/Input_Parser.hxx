#pragma once

#include "Polynomial_Vector_Matrix_State.hxx"

using namespace std::string_literals;
class Input_Parser : public xmlpp::SaxParser
{
public:
  std::string sdp_name="sdp";
  bool inside_sdp = false, finished_objective = false,
       finished_polynomial_vector_matrices = false;

  Vector_State<Element_State> objective_state;
  Vector_State<Polynomial_Vector_Matrix_State> polynomial_vector_matrices_state;

  Input_Parser()
      : objective_state({"objective"s, "elt"s}),
        polynomial_vector_matrices_state(
          {"polynomialVectorMatrices"s, "polynomialVectorMatrix"s})
  {}
  ~Input_Parser() override = default;

protected:
  void on_start_element(const Glib::ustring &name,
                        const AttributeList &properties) override;
  void on_end_element(const Glib::ustring &name) override;
  void on_characters(const Glib::ustring &characters) override;
};
