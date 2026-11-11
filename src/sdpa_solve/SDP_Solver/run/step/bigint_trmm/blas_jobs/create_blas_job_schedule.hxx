#pragma once

#include "Blas_Job.hxx"
#include "sdpb_util/bigint_shared_memory/blas_jobs/Blas_Job_Schedule.hxx"
#include "sdpb_util/Verbosity.hxx"

namespace Sdpb::Sdpa::Trmm
{
  // Jobs for X =: L X (or X := X L) modulo each prime
  // L is Block_Diagonal_Matrix
  // X is vector<Block_Diagonal_Matrix>, vector length = vector_size
  // L and each element of X have the same block structure defined by block_dims
  //
  // If side == El::LEFT
  //   X := L X
  //   All elements of vector X are stacked horizontally,
  //   i.e. each block of X is now wide, ~ M_i x (M_i*vector_size)
  // If side == El::RIGHT:
  //   X := X L
  //   All elements of vector X are stacked vertically,
  //   i.e. each block of X is now tall, ~ (M_i * vector_size) x M_i
  // Each block of L_i is triangular, ~ M_i x M_i.
  Blas_Job_Schedule<Trmm::Blas_Job>
  create_blas_job_schedule(El::LeftOrRight side, El::UpperOrLower uplo,
                           El::Orientation orientation, El::UnitOrNonUnit diag,
                           size_t num_ranks, size_t num_primes,
                           const std::vector<size_t> &heights,
                           const std::vector<size_t> &widths,
                           Verbosity verbosity);
}