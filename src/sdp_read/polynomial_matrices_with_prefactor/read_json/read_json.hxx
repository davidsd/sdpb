#pragma once

#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <El.hpp>

#include <filesystem>
#include <vector>

void read_json(
  const std::filesystem::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix,
  std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization, size_t &num_matrices,
  std::map<size_t, Positive_Matrix_With_Prefactor> &parsed_matrices);