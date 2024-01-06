#pragma once

#include "polynomial_matrices_with_prefactor/PMWP_SDP.hxx"
#include "read_pvm_input.hxx"

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
