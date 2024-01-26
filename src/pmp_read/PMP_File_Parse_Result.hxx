#pragma once

#include "pmp/Polynomial_Vector_Matrix.hxx"

#include <El.hpp>

#include <filesystem>
#include <vector>

struct PMP_File_Parse_Result
{
  std::filesystem::path path;
  // Vector a_0..a_N, see (3.1) in SDPB Manual
  // Empty if no objectives in file
  std::vector<El::BigFloat> objective;
  // Normaliation vector n_0..n_N, see (3.1) in SDPB Manual
  // Empty if no normalization in file
  std::vector<El::BigFloat> normalization;
  // Total number of PMWP matrices in file
  size_t num_matrices = 0;
  // If file is read by several processes,
  // each process saves only some matrices, according to should_parse_matrix()
  // parsed_matrices is a map: index -> matrix,
  // where 0 <= index < num_matrices
  std::map<size_t, Polynomial_Vector_Matrix> parsed_matrices;

  PMP_File_Parse_Result() = default;

  static void validate(const PMP_File_Parse_Result &result);

  static PMP_File_Parse_Result
  read(const std::filesystem::path &input_path, bool should_parse_objective,
       bool should_parse_normalization,
       const std::function<bool(size_t matrix_index)> &should_parse_matrix);

  // Allow moving and prevent accidential copying

  PMP_File_Parse_Result(const PMP_File_Parse_Result &other) = delete;
  PMP_File_Parse_Result(PMP_File_Parse_Result &&other) noexcept = default;
  PMP_File_Parse_Result &operator=(const PMP_File_Parse_Result &other)
    = delete;
  PMP_File_Parse_Result &operator=(PMP_File_Parse_Result &&other) noexcept
    = default;
};
