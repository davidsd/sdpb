#pragma once

#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <El.hpp>

#include <boost/noncopyable.hpp>

#include <filesystem>
#include <vector>

// PMWP is Positive_Matrix_With_Prefactor
struct PMWP_File_Parse_Result : boost::noncopyable
{
  std::filesystem::path path;
  // Empty if no objectives in file
  std::vector<El::BigFloat> objective;
  // Empty if no normalization in file
  std::vector<El::BigFloat> normalization;
  // Total number of PMWP matrices in file
  size_t num_matrices;
  // If file is read by several processes,
  // each process saves only some matrices, according to should_parse_matrix()
  // parsed_matrices is a map: index -> matrix,
  // where 0 <= index < num_matrices
  std::map<size_t, Positive_Matrix_With_Prefactor> parsed_matrices;

  PMWP_File_Parse_Result(
    const std::filesystem::path &input_path,
    const std::function<bool(size_t matrix_index)> &should_parse_matrix);
};