#pragma once

#include "pmp2sdp/write_sdp.hxx"
#include "sdpb_util/Number_State.hxx"
#include "sdpb_util/Vector_State.hxx"

using namespace std::string_literals;
class Xml_Polynomial_Vector_Matrix_State
{
public:
  std::string name;
  std::string rows_name = "rows", columns_name = "cols", rows_string,
              columns_string;
  bool inside = false, inside_rows = false, inside_columns = false;
  int rows = 0, cols = 0;
  std::unique_ptr<Polynomial_Vector_Matrix> value;
  std::function<void(Xml_Polynomial_Vector_Matrix_State &matrix_state)>
    process_matrix;

  using Polynomial_State = Vector_State<Number_State<El::BigFloat>>;
  using Polynomial_Vector_State = Vector_State<Polynomial_State>;
  Vector_State<Polynomial_Vector_State> elements_state;
  Vector_State<Number_State<El::BigFloat>> sample_points_state;
  Vector_State<Number_State<El::BigFloat>> sample_scalings_state;
  Vector_State<Polynomial_State> bilinear_basis_state;

  Xml_Polynomial_Vector_Matrix_State(
    const std::vector<std::string> &names, const size_t &offset,
    const std::function<void(Xml_Polynomial_Vector_Matrix_State &matrix_state)>
      &process_matrix)
      : name(names.at(offset)),
        process_matrix(process_matrix),
        elements_state(
          {"elements"s, "polynomialVector"s, "polynomial"s, "coeff"s}),
        sample_points_state({"samplePoints"s, "elt"s}),
        sample_scalings_state({"sampleScalings"s, "elt"s}),
        bilinear_basis_state({"bilinearBasis"s, "polynomial"s, "coeff"s})
  {}

  bool xml_on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        if(element_name == rows_name)
          {
            inside_rows = true;
          }
        else if(element_name == columns_name)
          {
            inside_columns = true;
          }
        else if(!(elements_state.xml_on_start_element(element_name)
                  || sample_points_state.xml_on_start_element(element_name)
                  || sample_scalings_state.xml_on_start_element(element_name)
                  || bilinear_basis_state.xml_on_start_element(element_name)))
          {
            RUNTIME_ERROR(
              "Inside polynomialVectorMatrix, expected 'rows', 'cols', "
              "'elements', 'sample_points', 'sample_scalings', or "
              "'bilinear_basis', but found '",
              element_name, "'");
          }
      }
    else if(element_name == name)
      {
        inside = true;
        elements_state.value.clear();
        sample_points_state.value.clear();
        sample_scalings_state.value.clear();
        bilinear_basis_state.value.clear();
        rows_string.clear();
        columns_string.clear();
      }
    return inside;
  }

  bool xml_on_end_element(const std::string &element_name)
  {
    bool result(inside);
    if(inside)
      {
        if(element_name == name)
          {
            inside = false;
            process_matrix(*this);
          }
        else if(inside_rows && element_name == rows_name)
          {
            inside_rows = false;
            rows = std::stoi(rows_string);
          }
        else if(inside_columns && element_name == columns_name)
          {
            inside_columns = false;
            cols = std::stoi(columns_string);
          }
        else if(elements_state.xml_on_end_element(element_name)) {}
        else if(sample_points_state.xml_on_end_element(element_name)) {}
        else if(sample_scalings_state.xml_on_end_element(element_name)) {}
        else if(bilinear_basis_state.xml_on_end_element(element_name)) {}
      }
    return result;
  }

  bool xml_on_characters(const xmlChar *characters, int length)
  {
    if(inside)
      {
        if(inside_rows)
          {
            rows_string.append(reinterpret_cast<const char *>(characters),
                               length);
          }
        else if(inside_columns)
          {
            columns_string.append(reinterpret_cast<const char *>(characters),
                                  length);
          }
        else
          {
            elements_state.xml_on_characters(characters, length)
              || sample_points_state.xml_on_characters(characters, length)
              || sample_scalings_state.xml_on_characters(characters, length)
              || bilinear_basis_state.xml_on_characters(characters, length);
          }
      }
    return inside;
  }
};
