#pragma once

#include "Positive_Matrix_With_Prefactor.hxx"
#include "../sdp_convert.hxx"

#include <boost/filesystem.hpp>

void read_input(const boost::filesystem::path &input_file,
                std::vector<El::BigFloat> &objectives,
                std::vector<El::BigFloat> &normalization,
                std::vector<Positive_Matrix_With_Prefactor> &matrices);

void read_input(
  const std::vector<boost::filesystem::path> &input_files,
  std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices,
  size_t &num_processed);

