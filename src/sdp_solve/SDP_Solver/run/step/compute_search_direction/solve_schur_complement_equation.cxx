#include "sdp_solve/SDP_Solver/run/step/Q_Solve_Cache.hxx"
#include "sdp_solve/SDP_Solver.hxx"
#include "sdpb_util/copy_matrix.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <optional>

// Solve the Schur complement equation for dx, dy.
//
// - As inputs, dx and dy are the residues r_x and r_y on the
//   right-hand side of the Schur complement equation.
// - As outputs, dx and dy are overwritten with the solutions of the
//   Schur complement equation.
//
// The equation is solved using the block-decomposition described in
// the manual.
//
// If Q lives on a sub-grid of the ranks (--qSolveRanks, see Q_Grids.hxx), dy is summed over all ranks by one AllReduce, solved
// on Q's grid and broadcast from rank 0. Otherwise, as before, the
// contributions to dy are sent to the owners of Q's [MC,MR] vector
// (QueueUpdate / ProcessQueues).
//
// timers (may be nullptr): sub-timers of the stages below.
// q_solve_cache (may be nullptr; --qSolveCache): the same solve with Q,
// keeping the copies of Q's diagonal blocks for the step's other solves.
void solve_schur_complement_equation(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy,
  Timers *timers, Q_Solve_Cache *q_solve_cache)
{
  std::optional<Scoped_Timer> phase;
  auto begin = [&](const char *name) {
    phase.reset();
    if(timers != nullptr)
      phase.emplace(*timers, name);
  };
  // dx = schur_complement^{-1}.dx
  begin("L_solve");
  lower_triangular_solve(schur_complement_cholesky, dx);

  begin("dy_sum");
  El::Matrix<El::BigFloat> dy_sum;
  Zeros(dy_sum, Q.Height(), 1);
  for(size_t block = 0; block < schur_off_diagonal.blocks.size(); ++block)
    {
      Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
           schur_off_diagonal.blocks[block], dx.blocks[block],
           El::BigFloat(1), dy.blocks[block]);

      // Locally sum contributions to dy
      for(int64_t row = 0; row < dy.blocks[block].LocalHeight(); ++row)
        {
          int64_t global_row(dy.blocks[block].GlobalRow(row));
          for(int64_t column = 0; column < dy.blocks[block].LocalWidth();
              ++column)
            {
              int64_t global_column(dy.blocks[block].GlobalCol(column));
              dy_sum(global_row, global_column)
                += dy.blocks[block].GetLocal(row, column);
            }
        }
    }

  // dy_dist = Q^{-1}.dy_dist
  auto solve_Q = [&](El::DistMatrix<El::BigFloat> &dy_dist) {
    begin("Q_solve");
    if(q_solve_cache != nullptr)
      q_solve_cache->solve(Q, dy_dist, timers);
    else
      El::cholesky::SolveAfter(El::UpperOrLowerNS::UPPER,
                               El::OrientationNS::NORMAL, Q, dy_dist);
  };
  // dx += schur_off_diagonal.dy
  auto update_dx = [&](const auto &dy_all) {
    begin("dx_update");
    for(size_t block = 0; block < schur_off_diagonal.blocks.size(); ++block)
      {
        copy_matrix(dy_all, dy.blocks[block]);
        Gemv(El::OrientationNS::NORMAL, El::BigFloat(1),
             schur_off_diagonal.blocks[block], dy.blocks[block],
             El::BigFloat(1), dx.blocks[block]);
      }
  };

  if(Q.Grid() != El::Grid::Default())
    {
      // Q on a sub-grid
      begin("dy_allreduce");
      El::mpi::AllReduce(dy_sum.Buffer(), dy_sum.Height(), El::mpi::SUM,
                         El::mpi::COMM_WORLD);
      if(Q.Participating())
        {
          El::DistMatrix<El::BigFloat> dy_dist(Q.Grid());
          dy_dist.Resize(Q.Height(), 1);
          copy_matrix(dy_sum, dy_dist);
          solve_Q(dy_dist);
          begin("dy_star_star");
          El::DistMatrix<El::BigFloat, El::STAR, El::STAR> dy_star_star(
            dy_dist);
          if(El::mpi::Rank() == 0)
            dy_sum = dy_star_star.Matrix();
        }
      begin("dy_bcast");
      El::mpi::Broadcast(dy_sum.Buffer(), dy_sum.Height(), 0,
                         El::mpi::COMM_WORLD);
      update_dx(dy_sum);
    }
  else
    {
      // Send out updates for dy
      begin("dy_queue");
      El::DistMatrix<El::BigFloat> dy_dist;
      Zeros(dy_dist, Q.Height(), 1);
      El::BigFloat zero(0);
      for(int64_t row = 0; row < dy_sum.Height(); ++row)
        for(int64_t column = 0; column < dy_sum.Width(); ++column)
          {
            if(dy_sum(row, column) != zero)
              {
                dy_dist.QueueUpdate(row, column, dy_sum(row, column));
              }
          }
      begin("dy_process_queues");
      dy_dist.ProcessQueues();

      solve_Q(dy_dist);
      begin("dy_star_star");
      El::DistMatrix<El::BigFloat, El::STAR, El::STAR> dy_local(dy_dist);
      update_dx(dy_local);
    }

  // dx = schur_complement^{-T}.dx
  begin("LT_solve");
  lower_triangular_transpose_solve(schur_complement_cholesky, dx);
  phase.reset();
}
