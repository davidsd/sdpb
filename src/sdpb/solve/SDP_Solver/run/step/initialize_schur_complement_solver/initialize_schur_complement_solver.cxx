#include "../../../../SDP.hxx"
#include "../../../../Block_Diagonal_Matrix.hxx"
#include "../../../../../../Timers.hxx"

// Compute the quantities needed to solve the Schur complement
// equation
//
// {{S, -B}, {B^T, 0}} . {dx, dy} = {r, s}
//
// using the method described in the manual:
//
// - S = Tr(A_X_Inv A_Y)
//
// - Cholesky decomposition S = L L^T.
//
// - L_inv_B = L^{-1} B
//
// - Q = L_inv_B^T L_inv_B
//
// - Cholesky decomposition of Q.

void reduce_and_scatter(const El::mpi::Comm &comm,
                        std::vector<El::byte> &send_buffer,
                        std::vector<El::byte> &receive_buffer,
                        std::vector<int> &rank_sizes);

void compute_S(const Block_Info &block_info,
               const Block_Diagonal_Matrix &A_X_inv,
               const Block_Diagonal_Matrix &A_Y, Block_Diagonal_Matrix &S,
               Timers &timers);

void initialize_Q_group(const SDP &sdp, const Block_Info &block_info,
                        const Block_Diagonal_Matrix &S, Block_Matrix &L_inv_B,
                        Block_Diagonal_Matrix &L,
                        El::DistMatrix<El::BigFloat> &Q_group, Timers &timers);

void fill_send_buffer(const El::DistMatrix<El::BigFloat> &Q,
                      const El::DistMatrix<El::BigFloat> &Q_group,
                      std::vector<El::byte> &send_buffer,
                      std::vector<int> &rank_sizes, size_t &serialized_size,
                      Timers &timers);

void synchronize_Q(std::vector<El::byte> &sending_buffer,
                   std::vector<int> &rank_sizes, const size_t &serialized_size,
                   El::DistMatrix<El::BigFloat> &Q, Timers &timers);

void initialize_schur_complement_solver(
  const Block_Info &block_info, const SDP &sdp,
  const Block_Diagonal_Matrix &A_X_inv, const Block_Diagonal_Matrix &A_Y,
  const El::Grid &group_grid, Block_Diagonal_Matrix &L, Block_Matrix &L_inv_B,
  El::DistMatrix<El::BigFloat> &Q, Timers &timers)
{
  auto &initialize_timer(
    timers.add_and_start("run.step.initializeSchurComplementSolver"));
  Block_Diagonal_Matrix S(block_info.schur_block_sizes,
                          block_info.block_indices,
                          block_info.schur_block_sizes.size(), group_grid);
  compute_S(block_info, A_X_inv, A_Y, S, timers);

  auto &Q_computation_timer(
    timers.add_and_start("run.step.initializeSchurComplementSolver.Q"));

  std::vector<El::byte> send_buffer;
  std::vector<int> rank_sizes(El::mpi::Size(El::mpi::COMM_WORLD));
  size_t serialized_size;
  {
    El::DistMatrix<El::BigFloat> Q_group(Q.Height(), Q.Width(), group_grid);
    initialize_Q_group(sdp, block_info, S, L_inv_B, L, Q_group, timers);
    fill_send_buffer(Q, Q_group, send_buffer, rank_sizes, serialized_size,
                     timers);
  }
  synchronize_Q(send_buffer, rank_sizes, serialized_size, Q, timers);
  Q_computation_timer.stop();

  auto &Cholesky_timer(
    timers.add_and_start("run.step.initializeSchurComplementSolver."
                         "Cholesky"));
  Cholesky(El::UpperOrLowerNS::UPPER, Q);
  Cholesky_timer.stop();
  initialize_timer.stop();
}
