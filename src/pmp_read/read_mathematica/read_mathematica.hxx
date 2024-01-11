#pragma once

#include "pmp/Polynomial_Vector_Matrix.hxx"

#include <El.hpp>

#include <filesystem>
#include <vector>

void read_mathematica(
  const std::filesystem::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix,
  std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization, size_t &num_matrices,
  std::map<size_t, Polynomial_Vector_Matrix> &parsed_matrices);