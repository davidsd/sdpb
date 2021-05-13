#include "../../sdp_solve.hxx"

// TODO: Have this be part of sdp_solve.hxx
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
                   const Block_Vector &y, const Block_Diagonal_Matrix &Y,
                   Block_Vector &dx, Block_Vector &dy)
{
  std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    Q_Y_Q;

  compute_Q_Y_Q(block_info, Y, d_sdp.bases_blocks, Q_Y_Q);

  const El::BigFloat mu("2.620639431472500199491970787109003213463978122955701938670174937806705995823423358021528042236620972427419161662933951265307808522014595283650526821323883163144965569169316749861052114026692552518230484184675643194883464709791913895161981921897887713631347326476983452838785936273639286638657793231503534625961e-31");
  const El::BigFloat inv_sqrt_mu(1/El::Sqrt(mu));
  
  // const El::BigFloat inv_sqrt_mu(1e10);
  for(auto &parity: Q_Y_Q)
    for(auto &block: parity)
      for(auto &row: block)
        for(auto &column: row)
          {
            El::Scale(inv_sqrt_mu, column);
          }
  
  for(size_t block_index(0); block_index != dx.blocks.size(); ++block_index)
    {
      // dx = (c + dc) - (B + dB).y
      dx.blocks[block_index] = d_sdp.primal_objective_c.blocks[block_index];
      // dx.blocks[block_index] = sdp.primal_objective_c.blocks[block_index];
      // dx.blocks[block_index] += d_sdp.primal_objective_c.blocks[block_index];

      El::DistMatrix BdB(d_sdp.free_var_matrix.blocks[block_index]);
      // El::DistMatrix BdB(sdp.free_var_matrix.blocks[block_index]);
      // BdB+=d_sdp.free_var_matrix.blocks[block_index];
      Gemv(El::Orientation::NORMAL, El::BigFloat(-1),
           BdB, y.blocks[block_index],
           El::BigFloat(1), dx.blocks[block_index]);
      
      // dy = db - x.dB
      El::Zero(dy.blocks[block_index]);
      if(block_info.block_indices[block_index] == 0)
        {
          dy.blocks[block_index] = d_sdp.dual_objective_b;
        }
      El::Gemv(El::OrientationNS::TRANSPOSE, El::BigFloat(-1),
               d_sdp.free_var_matrix.blocks[block_index], x.blocks[block_index], El::BigFloat(1.0),
               dy.blocks[block_index]);
    }

  Block_Diagonal_Matrix schur_complement_cholesky(
    block_info.schur_block_sizes(), block_info.block_indices,
    block_info.num_points.size(), grid);
  Block_Matrix schur_off_diagonal;
  El::DistMatrix<El::BigFloat> Q(d_sdp.dual_objective_b.Height(),
                                 d_sdp.dual_objective_b.Height());

  Timers timers(false);
  // Repurpose this to compute A_Y_A_Y instead of A_X^-1_A_Y
  initialize_schur_complement_solver(block_info, sdp, Q_Y_Q, Q_Y_Q, grid,
                                     schur_complement_cholesky,
                                     schur_off_diagonal, Q, timers);

  // Solve for dx, dy in-place
  solve_schur_complement_equation(schur_complement_cholesky,
                                  schur_off_diagonal, Q, dx, dy);
  El::Print(dy.blocks.at(0),"dy");
  std::cout << "\n";
  El::Print(dx.blocks.at(0),"dx");
  std::cout << "\n";
}
