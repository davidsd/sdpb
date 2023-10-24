#include "sdp_solve/SDP.hxx"
#include "sdp_solve/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/Timers/Timers.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/Matrix_Normalizer.hxx"

// schur_off_diagonal = L^{-1} B
void initialize_schur_off_diagonal(
  const SDP &sdp, const Block_Info &block_info,
  const Block_Diagonal_Matrix &schur_complement,
  Block_Matrix &schur_off_diagonal,
  Block_Diagonal_Matrix &schur_complement_cholesky, Timers &timers)
{
  schur_off_diagonal.blocks.clear();
  schur_off_diagonal.blocks.reserve(schur_complement_cholesky.blocks.size());

  size_t num_blocks_local = schur_complement_cholesky.blocks.size();
  for(size_t block = 0; block < num_blocks_local; ++block)
    {
      auto block_index_string
        = std::to_string(block_info.block_indices[block]);
      {
        Scoped_Timer cholesky_timer(timers, "cholesky_" + block_index_string);
        schur_complement_cholesky.blocks[block]
          = schur_complement.blocks[block];

        Cholesky(El::UpperOrLowerNS::LOWER,
                 schur_complement_cholesky.blocks[block]);
      }

      // schur_off_diagonal = L^{-1} B
      Scoped_Timer solve_timer(timers, "solve_" + block_index_string);

      schur_off_diagonal.blocks.push_back(sdp.free_var_matrix.blocks[block]);
      El::Trsm(El::LeftOrRightNS::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
               El::BigFloat(1), schur_complement_cholesky.blocks[block],
               schur_off_diagonal.blocks[block]);
    }
}

// Q = P^T P = (L^{-1} B)^T (L^{-1} B) = schur_off_diagonal^T schur_off_diagonal
void initialize_Q_group(const SDP &sdp, const Block_Info &block_info,
                        Block_Matrix &schur_off_diagonal,
                        BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
                        El::DistMatrix<El::BigFloat> &Q_group, Timers &timers)
{
  // Explicitly deallocate the lower half of Q_group.  This
  // significantly reduces the total amount of memory required.
  El::Matrix<El::BigFloat> &local(Q_group.Matrix());
  for(int64_t row = 0; row < Q_group.Height(); ++row)
    for(int64_t column = 0; column < row; ++column)
      {
        if(Q_group.IsLocal(row, column))
          {
            mpf_clear(local(Q_group.LocalRow(row), Q_group.LocalCol(column))
                        .gmp_float.get_mpf_t());
            local(Q_group.LocalRow(row), Q_group.LocalCol(column))
              .gmp_float.get_mpf_t()[0]
              ._mp_d
              = nullptr;
          }
      }

  int block_width = Q_group.Width();

  std::vector<El::DistMatrix<El::BigFloat>> &P_blocks
    = schur_off_diagonal.blocks;
  Scoped_Timer normalizer_ctor_timer(timers, "Matrix_Normalizer_ctor");
  Matrix_Normalizer normalizer(P_blocks, block_width, El::gmp::Precision(),
                               bigint_syrk_context.shared_memory_comm);
  normalizer_ctor_timer.stop();

  {
    Scoped_Timer normalizer_timer(timers, "normalize_P");
    normalizer.normalize_and_shift_P_blocks(P_blocks);
  }
  auto uplo = El::UPPER;
  bigint_syrk_context.bigint_syrk_blas(uplo, P_blocks, Q_group, timers);

  Scoped_Timer restore_timer(timers, "restore_P_Q");
  normalizer.restore_P_blocks(P_blocks);

  // TODO: synchronize_Q, check that Q(i,i)==2^2N, then restore
  normalizer.restore_Q(uplo, Q_group);
}

void synchronize_Q(El::DistMatrix<El::BigFloat> &Q,
                   const El::DistMatrix<El::BigFloat> &Q_group,
                   Timers &timers);

void compute_Q(const SDP &sdp, const Block_Info &block_info,
               const Block_Diagonal_Matrix &schur_complement,
               Block_Matrix &schur_off_diagonal,
               Block_Diagonal_Matrix &schur_complement_cholesky,
               BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
               El::DistMatrix<El::BigFloat> &Q, Timers &timers)
{
  Scoped_Timer timer(timers, "Q");

  initialize_schur_off_diagonal(sdp, block_info, schur_complement,
                                schur_off_diagonal, schur_complement_cholesky,
                                timers);

  const El::Grid grid(bigint_syrk_context.shared_memory_comm);
  El::DistMatrix<El::BigFloat> Q_group(Q.Height(), Q.Width(), grid);
  initialize_Q_group(sdp, block_info, schur_off_diagonal, bigint_syrk_context,
                     Q_group, timers);
  synchronize_Q(Q, Q_group, timers);
}