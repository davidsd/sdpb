#pragma once

#include "../sdp_convert.hxx"

#include <boost/filesystem.hpp>

void read_pvm_input(
  const std::vector<boost::filesystem::path> &input_files,
  std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices,
  size_t &num_processed);

