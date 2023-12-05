#pragma once

#include "sdp_convert/sdp_convert.hxx"
#include "sdpb_util/Number_State.hxx"
#include "sdpb_util/Vector_State.hxx"

using namespace std::string_literals;
class Polynomial_Vector_Matrix_State
{
public:
  std::string name;
  std::string rows_name = "rows", columns_name = "cols", rows_string,
              columns_string;
  bool inside = false, inside_rows = false, inside_columns = false;
  Polynomial_Vector_Matrix value;
  std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices;
  const size_t rank = El::mpi::Rank(),
               num_procs = El::mpi::Size(El::mpi::COMM_WORLD);
  size_t &num_processed;

  using Polynomial_State = Vector_State<Number_State<El::BigFloat>>;
  using Polynomial_Vector_State = Vector_State<Polynomial_State>;
  Vector_State<Polynomial_Vector_State> elements_state;
  Vector_State<Number_State<El::BigFloat>> sample_points_state;
  Vector_State<Number_State<El::BigFloat>> sample_scalings_state;
  Vector_State<Polynomial_State> bilinear_basis_state;

  Polynomial_Vector_Matrix_State(
    const std::vector<std::string> &names, const size_t &offset,
    std::vector<Polynomial_Vector_Matrix> &Polynomial_vector_matrices,
    size_t &Num_processed)
      : name(names.at(offset)), polynomial_vector_matrices(Polynomial_vector_matrices),
        num_processed(Num_processed),
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
            // Jump through hoops so that we immediately clear the
            // polynomial_vector_matrix after constructing any needed
            // dual_constraint_groups.  This significantly reduces the
            // memory usage, but does complicate the code.
            if(num_processed % num_procs == rank)
              {
                polynomial_vector_matrices.emplace_back(value);
              }
            ++num_processed;
            value.clear();
          }
        else if(inside_rows && element_name == rows_name)
          {
            inside_rows = false;
            value.rows = std::stoi(rows_string);
          }
        else if(inside_columns && element_name == columns_name)
          {
            inside_columns = false;
            value.cols = std::stoi(columns_string);
          }
        else if(elements_state.xml_on_end_element(element_name))
          {
            if(!elements_state.inside)
              {
                swap(value.elements, elements_state.value);
              }
          }
        else if(sample_points_state.xml_on_end_element(element_name))
          {
            if(!sample_points_state.inside)
              {
                std::swap(value.sample_points, sample_points_state.value);
              }
          }
        else if(sample_scalings_state.xml_on_end_element(element_name))
          {
            if(!sample_scalings_state.inside)
              {
                std::swap(value.sample_scalings, sample_scalings_state.value);
              }
          }
        else if(bilinear_basis_state.xml_on_end_element(element_name))
          {
            if(!bilinear_basis_state.inside)
              {
                swap(value.bilinear_basis, bilinear_basis_state.value);
              }
          }
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
