#pragma once

#include <vector>
#include "Blas_Job.hxx"
#include "Blas_Job_Schedule.hxx"

// Create BLAS jobs for matrix multiplication Q = P^T P.
// Each job calculates a submatrix Q_IJ = P_I^T P_J modulo some prime
// Since Q is symmetric, jobs are created only for the upper half (I<=J).
// Jobs are distributed across num_ranks ranks on a node.
// Q is a NxN matrix, N=output_matrix_height
Blas_Job_Schedule
create_blas_job_schedule(size_t num_ranks, size_t num_primes,
                         El::Int output_matrix_height, bool debug);
