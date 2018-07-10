#pragma once

#include "../../Polynomial_Vector_Matrix.hxx"
#include "Element_State.hxx"
#include "Vector_State.hxx"

namespace
{
  // Convenience functions to avoid copies
  void swap(std::vector<Polynomial> &polynomials,
            std::vector<std::vector<El::BigFloat>> &elements_vector)
  {
    for(auto &elements : elements_vector)
      {
        polynomials.emplace_back();
        std::swap(polynomials.back().coefficients, elements);
      }
  }

  void swap(std::vector<std::vector<Polynomial>> &polynomials_vector,
            std::vector<std::vector<std::vector<El::BigFloat>>>
              &elements_vector_vector)
  {
    for(auto &elements_vector : elements_vector_vector)
      {
        polynomials_vector.emplace_back();
        swap(polynomials_vector.back(), elements_vector);
      }
  }
}

using namespace std::string_literals;
class Polynomial_Vector_Matrix_State
{
public:
  std::string name;
  bool inside = false, inside_rows = false, inside_columns = false;
  Polynomial_Vector_Matrix value;

  using Polynomial_State = Vector_State<Element_State>;
  using Polynomial_Vector_State = Vector_State<Polynomial_State>;
  Vector_State<Polynomial_Vector_State> elements_state;
  Vector_State<Element_State> sample_points_state;
  Vector_State<Element_State> sample_scalings_state;
  Vector_State<Polynomial_State> bilinear_basis_state;

  Polynomial_Vector_Matrix_State(const std::vector<std::string> &names,
                                 const size_t &offset)
      : name(names.at(offset)),
        elements_state(
          {"elements"s, "polynomialVector"s, "polynomial"s, "coeff"s}),
        sample_points_state({"samplePoints"s, "elt"s}),
        sample_scalings_state({"sampleScalings"s, "elt"s}),
        bilinear_basis_state({"bilinearBasis"s, "polynomial"s, "coeff"s})
  {}

  bool on_start_element(const Glib::ustring &element_name,
                        const xmlpp::SaxParser::AttributeList &attributes)
  {
    if(inside)
      {
        if(element_name == "rows")
          {
            inside_rows = true;
          }
        else if(element_name == "cols")
          {
            inside_columns = true;
          }
        else if(!(elements_state.on_start_element(element_name, attributes)
                  || sample_points_state.on_start_element(element_name,
                                                          attributes)
                  || sample_scalings_state.on_start_element(element_name,
                                                            attributes)
                  || bilinear_basis_state.on_start_element(element_name,
                                                           attributes)))
          {
            throw std::runtime_error(
              "Inside polynomialVectorMatrix, expected 'rows', 'cols', "
              "'elements', 'sample_points', 'sample_scalings', or "
              "'bilinear_basis', but found '"
              + element_name + "'");
          }
      }
    else if(element_name == name)
      {
        inside = true;
        elements_state.value.clear();
        sample_points_state.value.clear();
        sample_scalings_state.value.clear();
        bilinear_basis_state.value.clear();
      }
    return inside;
  }

  bool on_end_element(const Glib::ustring &element_name)
  {
    bool result(false);

    if(inside)
      {
        if(element_name == name)
          {
            inside = false;
            result = true;
          }
        else if(element_name == "rows")
          {
            inside_rows = false;
          }
        else if(element_name == "cols")
          {
            inside_columns = false;
          }
        else if(elements_state.on_end_element(element_name))
          {
            swap(value.elements, elements_state.value);
          }
        else if(sample_points_state.on_end_element(element_name))
          {
            std::swap(value.sample_points, sample_points_state.value);
          }
        else if(sample_scalings_state.on_end_element(element_name))
          {
            std::swap(value.sample_scalings, sample_scalings_state.value);
          }
        else if(bilinear_basis_state.on_end_element(element_name))
          {
            swap(value.bilinear_basis, bilinear_basis_state.value);
          }
      }
    return result;
  }

  bool on_characters(const Glib::ustring &characters)
  {
    if(inside)
      {
        if(inside_rows)
          {
            value.rows = std::stoi(characters);
          }
        else if(inside_columns)
          {
            value.cols = std::stoi(characters);
          }
        else
          {
            elements_state.on_characters(characters)
              || sample_points_state.on_characters(characters)
              || sample_scalings_state.on_characters(characters)
              || bilinear_basis_state.on_characters(characters);
          }
      }
    return inside;
  }
};
