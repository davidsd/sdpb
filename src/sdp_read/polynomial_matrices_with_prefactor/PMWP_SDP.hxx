#pragma once

#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <El.hpp>
#include <boost/noncopyable.hpp>

#include <filesystem>

// PMWP is Positive_Matrix_With_Prefactor
// Format described in SDPB Manual, eq. (3.1)
struct PMWP_SDP : boost::noncopyable
{
  // vector a_i, i=0..N
  std::vector<El::BigFloat> objective;
  // normalization vector n_i, i=0..N
  std::vector<El::BigFloat> normalization;
  // Total number of PMWP matrices
  size_t num_matrices;
  // In case of several processes,
  // each process owns only some matrices.
  std::vector<Positive_Matrix_With_Prefactor> matrices;
  // global index of matrices[i], lies in [0..num_matrices)
  std::vector<size_t> matrix_index_local_to_global;

  explicit PMWP_SDP(const std::filesystem::path &input_path);
  explicit PMWP_SDP(const std::vector<std::filesystem::path> &input_paths);
};
