#pragma once

#include "sdp_read/read_input.hxx"
#include "sdp_read/read_pvm_input.hxx"

std::vector<boost::filesystem::path>
read_file_list(const boost::filesystem::path &input_file);

std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational);

std::vector<Boost_Float> sample_points(const size_t &num_points);
