#include "../../../../SDP.hxx"
#include "../../../../Block_Diagonal_Matrix.hxx"
#include "../../../../../../Timers.hxx"

// Compute the quantities needed to solve the Schur complement
// equation
//
// {{S, -B}, {B^T, 0}} . {dx, dy} = {r, s}
//
// (where S = SchurComplement, B = FreeVarMatrix), using the method
// described in the manual:
//
// - Compute S using BilinearPairingsXInv and BilinearPairingsY.
//
// - Compute the Cholesky decomposition S' = L' L'^T.
//
// - Form B' = (B U) and compute
//
//   - SchurOffDiagonal = L'^{-1} B
//   - L'^{-1} U
//   - Q = (L'^{-1} B')^T (L'^{-1} B') - {{0, 0}, {0, 1}}
//
// - Compute the LU decomposition of Q.
//
// This data is sufficient to efficiently solve the above equation for
// a given r,s.
//
// Inputs:
// - BilinearPairingsXInv, BilinearPairingsY (these are members of
//   SDPSolver, but we include them as arguments to emphasize that
//   they must be computed first)
// Workspace (members of SDPSolver which are modified by this method
// and not used later):
// - SchurComplement
// Outputs (members of SDPSolver which are modified by this method and
// used later):
// - SchurComplementCholesky
// - SchurOffDiagonal
//

void reduce_and_scatter(const El::mpi::Comm &comm,
                        std::vector<El::byte> &send_buffer,
                        std::vector<El::byte> &receive_buffer,
                        std::vector<int> &rank_sizes);

void compute_schur_complement(
  const Block_Info &block_info,
  const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y,
  Block_Diagonal_Matrix &schur_complement, Timers &timers);

void initialize_Q_group(const SDP &sdp, const Block_Info &block_info,
                        const Block_Diagonal_Matrix &schur_complement,
                        Block_Matrix &schur_off_diagonal,
                        Block_Diagonal_Matrix &schur_complement_cholesky,
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
  const Block_Diagonal_Matrix &bilinear_pairings_X_inv,
  const Block_Diagonal_Matrix &bilinear_pairings_Y, const El::Grid &group_grid,
  Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal, El::DistMatrix<El::BigFloat> &Q,
  Timers &timers)
{
  auto &initialize_timer(
    timers.add_and_start("run.step.initializeSchurComplementSolver"));
  // The Schur complement matrix S: a Block_Diagonal_Matrix with one
  // block for each 0 <= j < J.  SchurComplement.blocks[j] has dimension
  // (d_j+1)*m_j*(m_j+1)/2
  //
  Block_Diagonal_Matrix schur_complement(
    block_info.schur_block_sizes, block_info.block_indices,
    block_info.schur_block_sizes.size(), group_grid);

  compute_schur_complement(block_info, bilinear_pairings_X_inv,
                           bilinear_pairings_Y, schur_complement, timers);

  auto &Q_computation_timer(
    timers.add_and_start("run.step.initializeSchurComplementSolver.Q"));

  std::vector<El::byte> send_buffer;
  std::vector<int> rank_sizes(El::mpi::Size(El::mpi::COMM_WORLD));
  size_t serialized_size;
  {
    El::DistMatrix<El::BigFloat> Q_group(Q.Height(), Q.Width(), group_grid);
    initialize_Q_group(sdp, block_info, schur_complement, schur_off_diagonal,
                       schur_complement_cholesky, Q_group, timers);
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
