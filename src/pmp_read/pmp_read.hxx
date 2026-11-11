#pragma once
#include "pmp/Polynomial_Matrix_Program.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>
#include <functional>
#include <optional>
#include <vector>

// Input files for which known_num_matrices(path) returns a value are not
// parsed: their matrices only take up block indices (used by pmp2sdp
// --reuseSdp). If file_matrix_counts is not null, it receives every input
// file with its number of matrices, in block order.
Polynomial_Matrix_Program read_polynomial_matrix_program(
  const Environment &env,
  const std::vector<std::filesystem::path> &input_files, int64_t max_num_poles,
  const Verbosity &verbosity, Timers &timers,
  const std::function<std::optional<size_t>(const std::filesystem::path &)>
    &known_num_matrices = nullptr,
  std::vector<std::pair<std::filesystem::path, size_t>> *file_matrix_counts
  = nullptr);

Polynomial_Matrix_Program read_polynomial_matrix_program(
  const Environment &env, const std::filesystem::path &input_file,
  int64_t max_num_poles, const Verbosity &verbosity, Timers &timers,
  const std::function<std::optional<size_t>(const std::filesystem::path &)>
    &known_num_matrices = nullptr,
  std::vector<std::pair<std::filesystem::path, size_t>> *file_matrix_counts
  = nullptr);

std::vector<std::filesystem::path>
read_nsv_file_list(const std::filesystem::path &input_file);

std::vector<std::filesystem::path>
collect_files_expanding_nsv(const std::filesystem::path &input_file);

std::vector<std::filesystem::path> collect_files_expanding_nsv(
  const std::vector<std::filesystem::path> &input_files);
