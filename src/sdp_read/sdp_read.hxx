#pragma once

#include "read_input.hxx"
#include "read_pvm_input.hxx"

std::vector<std::filesystem::path>
read_nsv_file_list(const std::filesystem::path &input_file);

std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational);

std::vector<Boost_Float> sample_points(const size_t &num_points);
