#pragma once
#include "pmp/Polynomial_Matrix_Program.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>
#include <vector>

Polynomial_Matrix_Program read_polynomial_matrix_program(
  const Environment &env,
  const std::vector<std::filesystem::path> &input_files,
  const Verbosity &verbosity, Timers &timers);

Polynomial_Matrix_Program
read_polynomial_matrix_program(const Environment &env,
                               const std::filesystem::path &input_file,
                               const Verbosity &verbosity, Timers &timers);

std::vector<std::filesystem::path>
read_nsv_file_list(const std::filesystem::path &input_file);

std::vector<std::filesystem::path>
collect_files_expanding_nsv(const std::filesystem::path &input_file);

std::vector<std::filesystem::path> collect_files_expanding_nsv(
  const std::vector<std::filesystem::path> &input_files);
