#pragma once

#include "sdp_convert/Polynomial_Vector_Matrix.hxx"
#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <El.hpp>

#include <boost/noncopyable.hpp>

#include <filesystem>
#include <vector>

// PVM is Polynomial_Vector_Matrix
struct PVM_File_Parse_Result : boost::noncopyable
{
  std::filesystem::path path;
  // {b_0..b_N}
  // Empty if no objectives in file
  std::vector<El::BigFloat> objective;
  // Total number of PMWP matrices in file
  size_t num_matrices;
  // If file is read by several processes,
  // each process saves only some matrices, according to should_parse_matrix()
  // parsed_matrices is a map: index -> matrix,
  // where 0 <= index < num_matrices
  std::map<size_t, Polynomial_Vector_Matrix> parsed_matrices;

  PVM_File_Parse_Result(
    const std::filesystem::path &input_path,
    const std::function<bool(size_t matrix_index)> &should_parse_matrix);
};