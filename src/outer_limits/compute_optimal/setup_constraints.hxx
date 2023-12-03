#pragma once

#include "outer_limits/Function.hxx"

void setup_constraints(
  const size_t &max_index, const size_t &num_blocks,
  const El::BigFloat &epsilon, const El::BigFloat &infinity,
  const std::vector<std::vector<std::vector<std::vector<Function>>>>
    &function_blocks,
  const std::vector<El::BigFloat> &normalization,
  const std::vector<std::set<El::BigFloat>> &points,
  std::vector<std::vector<El::BigFloat>> &primal_objective_c,
  std::vector<El::Matrix<El::BigFloat>> &free_var_matrix);
