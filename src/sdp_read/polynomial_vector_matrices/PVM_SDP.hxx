#pragma once

#include "sdp_convert/Polynomial_Vector_Matrix.hxx"
#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <El.hpp>
#include <boost/noncopyable.hpp>

#include<filesystem>

// PVM_SDP is Polynomial_Vector_Matrix
// Format described in SDPB Manual, eq. (2.2)
struct PVM_SDP : boost::noncopyable
{
  // vector b_i, i=0..N
  // note that it includes constant b_0
  std::vector<El::BigFloat> objective;
  // Total number of PVM matrices
  size_t num_matrices;
  // In case of several processes,
  // each process owns only some matrices.
  std::vector<Polynomial_Vector_Matrix> matrices;
  // global index of matrices[i], lies in [0..num_matrices)
  std::vector<size_t> matrix_index_local_to_global;

  explicit PVM_SDP(const std::filesystem::path& input_path);
  explicit PVM_SDP(const std::vector<std::filesystem::path>& input_paths);
};
