#pragma once

#include "sdpb_util/Environment.hxx"

#include <El.hpp>

#include <memory>

// A process grid on a subset of the ranks for the solves with the factor of
// the global Q, the N x N matrix of the Schur complement solver
// (--qSolveRanks).
//
// Q is computed and factored on all ranks. The two triangular solves per
// search direction are latency-bound: on many nodes they get slower as ranks
// are added (each of El's Trsv blocks gathers its 128 x 128 diagonal block to
// every rank of the grid). With --qSolveRanks s, the factor is copied to a
// grid of s ranks for the solves (dy is summed over all ranks by one
// AllReduce and broadcast back after the solve).
//
// The ranks of a sub-grid are the first r in (node, rank on node) order, so
// s = k * (ranks per node) gives the ranks of k whole nodes. The grid is
// h x (s / h) with h the largest divisor of s not above sqrt(s). World rank
// 0 is always in it.
//
// The solve grid's shape and the AllReduce of dy change the order of the
// partial sums in the solves, so --qSolveRanks changes the iterates at
// working precision.
class Q_Grids
{
public:
  // solve_ranks <= 0 or >= number of ranks: no sub-grid.
  Q_Grids(const Environment &env, int64_t solve_ranks);

  // Grid for the solves with Q's factor; nullptr: the grid Q was computed on.
  [[nodiscard]] const El::Grid *solve_grid() const;

private:
  std::unique_ptr<El::Grid> solve_;
};

// B := the upper triangle of A, with B's strictly lower triangle zero.
// A and B are [MC,MR] on grids whose viewing communicator has the ranks of
// COMM_WORLD (B keeps its grid); all ranks call it. One AllToAll moves the
// N(N+1)/2 entries that Cholesky(UPPER) and the triangular solves read,
// instead of El::Copy's N^2.
void copy_upper_between_grids(const El::DistMatrix<El::BigFloat> &A,
                              El::DistMatrix<El::BigFloat> &B);
