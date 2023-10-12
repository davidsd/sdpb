#pragma once

#include <vector>
#include "Blas_Job.hxx"

// Create BLAS jobs for matrix multiplication Q = P^T P.
// Each job calculates a submatrix Q_IJ = P_I^T P_J modulo some prime
// Since Q is symmetric, jobs are created only for the upper half (I<=J).
// Jobs are distributed among num_ranks ranks, see Blas_Job_Schedule.
// Q is a NxN matrix, N=output_matrix_height
std::vector<Blas_Job> create_blas_jobs(size_t num_ranks, size_t num_primes,
                                       El::Int output_matrix_height);