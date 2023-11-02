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

// Check that Q_ii = 2^2N, where N = normalizer.precision.
// This follows from the fact that columns of P are normalized and multiplied by 2^N
void check_normalized_Q_diagonal(El::DistMatrix<El::BigFloat> &Q,
                                 const Matrix_Normalizer &normalizer,
                                 Timers &timers)
{
  Scoped_Timer timer(timers, "check_diagonal");
  for(int iLoc = 0; iLoc < Q.LocalHeight(); ++iLoc)
    for(int jLoc = 0; jLoc < Q.LocalWidth(); ++jLoc)
      {
        int i = Q.GlobalRow(iLoc);
        int j = Q.GlobalCol(jLoc);
        if(i == j)
          {
            auto value = Q.GetLocal(iLoc, jLoc); // should be 2^2N
            auto should_be_one = value >> 2 * normalizer.precision;

            auto diff = El::Abs(should_be_one - El::BigFloat(1));
            // diff should be equal to zero up to some precision.
            // We cannot control rounding errors exactly,
            // so we (conservatively) require that at least N/2 bits are correct.
            auto eps = El::BigFloat(1) >> normalizer.precision / 2;
            if(diff > eps)
              {
                El::RuntimeError(
                  "Normalized Q should have ones on diagonal. For i = ", i,
                  ": Q_ii = ", should_be_one, ", |Q_ii - 1| = ", diff,
                  ", eps = ", eps);
              }
          }
      }
}

// Q = P^T P = (L^{-1} B)^T (L^{-1} B) = schur_off_diagonal^T schur_off_diagonal
void syrk_Q(Block_Matrix &schur_off_diagonal,
            BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
            El::DistMatrix<El::BigFloat> &Q, Timers &timers)
{
  std::vector<El::DistMatrix<El::BigFloat>> &P_blocks
    = schur_off_diagonal.blocks;

  // Normalize P columns and multiply by 2^N
  int block_width = Q.Width();
  Scoped_Timer normalizer_ctor_timer(timers, "Matrix_Normalizer_ctor");
  Matrix_Normalizer normalizer(P_blocks, block_width, El::gmp::Precision(),
                               bigint_syrk_context.shared_memory_comm);
  normalizer_ctor_timer.stop();
  {
    Scoped_Timer normalizer_timer(timers, "normalize_P");
    normalizer.normalize_and_shift_P_blocks(P_blocks);
  }

  // Calculate Q = P^T P
  auto uplo = El::UPPER;
  bigint_syrk_context.bigint_syrk_blas(uplo, P_blocks, Q, timers);

  // Check that Q_ii = 2^2N
  check_normalized_Q_diagonal(Q, normalizer, timers);

  // Remove normalization
  Scoped_Timer restore_timer(timers, "restore_P_Q");
  normalizer.restore_P_blocks(P_blocks);
  normalizer.restore_Q(uplo, Q);
}

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
  syrk_Q(schur_off_diagonal, bigint_syrk_context, Q, timers);
}