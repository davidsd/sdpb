#pragma once
#include "pmp/Polynomial_Matrix_Program.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>
#include <vector>

Polynomial_Matrix_Program read_polynomial_matrix_program(
  const Environment &env,
  const std::vector<std::filesystem::path> &input_files, Timers &timers);

Polynomial_Matrix_Program
read_polynomial_matrix_program(const Environment &env,
                               const std::filesystem::path &input_file,
                               Timers &timers);

std::vector<std::filesystem::path>
read_nsv_file_list(const std::filesystem::path &input_file);

std::vector<std::filesystem::path>
collect_files_expanding_nsv(const std::filesystem::path &input_file);

std::vector<std::filesystem::path> collect_files_expanding_nsv(
  const std::vector<std::filesystem::path> &input_files);

std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational);

std::vector<Boost_Float> sample_points(const size_t &num_points);
