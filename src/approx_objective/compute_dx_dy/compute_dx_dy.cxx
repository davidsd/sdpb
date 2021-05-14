#include "../../sdp_solve.hxx"

// TODO: Have this be part of sdp_solve.hxx
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L);
void compute_Q_X_inv_Q(
  const Block_Info &block_info, const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &Q_X_inv_Q);

void compute_Q_Y_Q(
  const Block_Info &block_info, const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &Q_Y_Q);

void compute_dual_residues_and_error(
  const Block_Info &block_info, const SDP &sdp, const Block_Vector &y,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &Q_Y_Q,
  Block_Vector &dual_residues, El::BigFloat &dual_error, Timers &timers);

void initialize_schur_complement_solver(
  const Block_Info &block_info, const SDP &sdp,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &Q_X_inv_Q,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &Q_Y_Q,
  const El::Grid &block_grid, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal, El::DistMatrix<El::BigFloat> &Q,
  Timers &timers);

void solve_schur_complement_equation(
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx, Block_Vector &dy);

void compute_dx_dy(const Block_Info &block_info, const El::Grid &grid,
                   const SDP &sdp, const SDP &d_sdp, const Block_Vector &x,
                   const Block_Vector &y, const Block_Diagonal_Matrix &X,
                   const Block_Diagonal_Matrix &Y, Block_Vector &dx,
                   Block_Vector &dy)
{
  std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    Q_X_inv_Q, Q_Y_Q;

  Block_Diagonal_Matrix X_cholesky(X);
  cholesky_decomposition(X, X_cholesky);
  compute_Q_X_inv_Q(block_info, X_cholesky, d_sdp.bases_blocks, Q_X_inv_Q);
  compute_Q_Y_Q(block_info, Y, d_sdp.bases_blocks, Q_Y_Q);

  for(size_t block_index(0); block_index != dx.blocks.size(); ++block_index)
    {
      // dx = -dc - dB.y
      dx.blocks[block_index] = d_sdp.primal_objective_c.blocks[block_index];
      Gemv(El::Orientation::NORMAL, El::BigFloat(1),
           d_sdp.free_var_matrix.blocks[block_index], y.blocks[block_index],
           El::BigFloat(-1), dx.blocks[block_index]);

      // dy = db - x.dB
      El::Zero(dy.blocks[block_index]);
      if(block_info.block_indices[block_index] == 0)
        {
          dy.blocks[block_index] = d_sdp.dual_objective_b;
        }
      El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
               d_sdp.free_var_matrix.blocks[block_index],
               x.blocks[block_index], El::BigFloat(1.0),
               dy.blocks[block_index]);
    }

  Block_Diagonal_Matrix schur_complement_cholesky(
    block_info.schur_block_sizes(), block_info.block_indices,
    block_info.num_points.size(), grid);
  Block_Matrix schur_off_diagonal;
  El::DistMatrix<El::BigFloat> Q(d_sdp.dual_objective_b.Height(),
                                 d_sdp.dual_objective_b.Height());

  Timers timers(false);
  initialize_schur_complement_solver(block_info, sdp, Q_X_inv_Q,
                                     // Q_Y_Q,
                                     Q_Y_Q, grid, schur_complement_cholesky,
                                     schur_off_diagonal, Q, timers);

  // Solve for dx, dy in-place
  solve_schur_complement_equation(schur_complement_cholesky,
                                  schur_off_diagonal, Q, dx, dy);
}
