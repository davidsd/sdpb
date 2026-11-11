#pragma once

#include <El.hpp>

#include <vector>

struct Timers;

// The Q solve, El::cholesky::SolveAfter(UPPER, NORMAL, U, x) with
// one right-hand side (two El::Trsv: UT with ADJOINT, then UN), as El does it
// line for line, except that the [STAR,STAR] copy of each diagonal block of
// U (Blocksize() = 128 rows) is made once per factor and kept, instead of
// once per Trsv. SDPB makes four Trsv per iteration with the same factor
// (predictor and corrector), and each copy is an all-rank gather of the
// grid. Same arithmetic, so the same bits.
//
// One cache per factor: with --qSolveCache, SDP_Solver::step makes it after
// Cholesky_Q (the copies by Node_Diag_Exchange, one per node) and drops it
// at the end of the step. node_shared_copies = false: El's copies (tests).
class Q_Solve_Cache
{
public:
  explicit Q_Solve_Cache(bool node_shared_copies);
  // x := (U^T U)^{-1} x, as El::cholesky::SolveAfter(UPPER, NORMAL, U, x)
  void solve(const El::DistMatrix<El::BigFloat> &U,
             El::DistMatrix<El::BigFloat> &x, Timers *timers = nullptr);
  // Build the copies now (solve() builds them on first use)
  void build(const El::DistMatrix<El::BigFloat> &U, Timers *timers = nullptr);

private:
  bool node_shared;
  const El::BigFloat *U_buffer = nullptr;
  El::Int U_height = -1;
  std::vector<El::Matrix<El::BigFloat>> blocks;
};
