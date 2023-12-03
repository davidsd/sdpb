#pragma once

#include "Positive_Matrix_With_Prefactor.hxx"
#include "sdp_convert/sdp_convert.hxx"

#include <filesystem>

void read_input(const std::filesystem::path &input_file,
                std::vector<El::BigFloat> &objectives,
                std::vector<El::BigFloat> &normalization,
                std::vector<Positive_Matrix_With_Prefactor> &matrices,
                size_t &num_processed);

inline void read_input(const std::filesystem::path &input_file,
                       std::vector<El::BigFloat> &objectives,
                       std::vector<El::BigFloat> &normalization,
                       std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  size_t num_processed(0);
  read_input(input_file, objectives, normalization, matrices, num_processed);
}

void read_input(
  const std::vector<std::filesystem::path> &input_files,
  std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices,
  size_t &num_processed);
