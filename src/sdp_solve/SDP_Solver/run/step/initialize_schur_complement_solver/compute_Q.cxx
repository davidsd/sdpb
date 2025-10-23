#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/Matrix_Normalizer.hxx"
#include "sdp_solve/memory_estimates.hxx"
#include "sdpb_util/Timers/Timers.hxx"

// schur_off_diagonal = L^{-1} B
void initialize_schur_off_diagonal(
  const Environment &env, const SDP &sdp, const Block_Info &block_info,
  const Block_Diagonal_Matrix &schur_complement,
  Block_Matrix &schur_off_diagonal,
  Block_Diagonal_Matrix &schur_complement_cholesky, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms, const Verbosity verbosity)
{
  schur_off_diagonal.blocks.clear();
  // Setting correct width is needed for Block_Matrix_Normalizer.
  // TODO: refactor Block_Matrix to ensure correct block width for all blocks
  schur_off_diagonal.width = sdp.dual_objective_b.Height();
  schur_off_diagonal.blocks.reserve(schur_complement_cholesky.blocks.size());

  size_t num_blocks_local = schur_complement_cholesky.blocks.size();
  for(size_t block = 0; block < num_blocks_local; ++block)
    {
      const auto global_block_index = block_info.block_indices.at(block);
      auto block_index_string = std::to_string(global_block_index);
      {
        Scoped_Timer cholesky_timer(timers, "cholesky_" + block_index_string);
        schur_complement_cholesky.blocks[block]
          = schur_complement.blocks[block];

        try
          {
            Cholesky(El::UpperOrLowerNS::LOWER,
                     schur_complement_cholesky.blocks[block]);
          }
        catch(std::exception &e)
          {
            RUNTIME_ERROR(
              "Error when computing Cholesky decomposition of block_",
              global_block_index, ": ", e.what());
          }
        block_timings_ms(global_block_index, 0)
          += cholesky_timer.elapsed_milliseconds();
      }

      // schur_off_diagonal = L^{-1} B
      Scoped_Timer solve_timer(timers, "solve_" + block_index_string);

      schur_off_diagonal.blocks.push_back(sdp.free_var_matrix.blocks[block]);
      El::Trsm(El::LeftOrRightNS::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
               El::BigFloat(1), schur_complement_cholesky.blocks[block],
               schur_off_diagonal.blocks[block]);
      block_timings_ms(global_block_index, 0)
        += solve_timer.elapsed_milliseconds();
    }
  if(verbosity >= Verbosity::trace)
    {
      print_allocation_message_per_node(
        env, "schur_off_diagonal", get_allocated_bytes(schur_off_diagonal));
      // schur_complement_cholesky was already allocated in SDP_Solver::step()
    }
}

// Check that Q_ii = 2^2N, where N = normalizer.precision.
// This follows from the fact that columns of P are normalized and multiplied by 2^N
void check_normalized_Q_diagonal(const El::DistMatrix<El::BigFloat> &Q,
                                 const Block_Matrix_Normalizer<> &normalizer,
                                 Timers &timers)
{
  Scoped_Timer timer(timers, "check_diagonal");
  for(int iLoc = 0; iLoc < Q.LocalHeight(); ++iLoc)
    for(int jLoc = 0; jLoc < Q.LocalWidth(); ++jLoc)
      {
        const int i = Q.GlobalRow(iLoc);
        const int j = Q.GlobalCol(jLoc);
        if(i == j)
          {
            auto value = Q.GetLocal(iLoc, jLoc); // should be 2^2N
            auto should_be_one = value >> 2 * normalizer.bits;

            auto diff = El::Abs(should_be_one - El::BigFloat(1));
            // diff should be equal to zero up to some precision.
            // We cannot control rounding errors exactly,
            // so we (conservatively) require that at least N/2 bits are correct.
            auto eps = El::BigFloat(1) >> normalizer.bits / 2;
            ASSERT(diff < eps,
                   "Normalized Q should have ones on diagonal. For i = ", i,
                   ": Q_ii = ", should_be_one, ", |Q_ii - 1| = ", diff,
                   ", eps = ", eps);
          }
      }
}

// Q = P^T P = (L^{-1} B)^T (L^{-1} B) = schur_off_diagonal^T schur_off_diagonal
void syrk_Q(const Environment &env, Block_Matrix &schur_off_diagonal,
            BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
            El::DistMatrix<El::BigFloat> &Q, Timers &timers,
            El::Matrix<int32_t> &block_timings_ms, const Verbosity verbosity)
{
  Scoped_Timer syrk_timer(timers, "syrk");

  // Normalize P columns and multiply by 2^N
  Scoped_Timer normalizer_timer(timers, "normalize_P");
  const auto normalizer
    = normalize_and_shift<Matrix_Normalization_Kind::COLUMNS>(
      schur_off_diagonal, El::gmp::Precision(), El::mpi::COMM_WORLD);
  if(verbosity >= Verbosity::trace)
    {
      print_allocation_message_per_node(env, "Matrix_Normalizer",
                                        get_allocated_bytes(normalizer.norms));
    }
  normalizer_timer.stop();

  // Calculate Q = P^T P
  constexpr auto uplo = El::UPPER;
  bigint_syrk_context.bigint_syrk_blas(uplo, schur_off_diagonal.blocks, Q,
                                       timers, block_timings_ms);

  // Check that Q_ii = 2^2N
  check_normalized_Q_diagonal(Q, normalizer, timers);

  // Remove normalization
  Scoped_Timer restore_timer(timers, "restore_P_Q");
  normalizer.restore(schur_off_diagonal);
  restore_syrk_output(uplo, normalizer, Q);
}

void compute_Q(const Environment &env, const SDP &sdp,
               const Block_Info &block_info,
               const Block_Diagonal_Matrix &schur_complement,
               Block_Matrix &schur_off_diagonal,
               Block_Diagonal_Matrix &schur_complement_cholesky,
               BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
               El::DistMatrix<El::BigFloat> &Q, Timers &timers,
               El::Matrix<int32_t> &block_timings_ms,
               const Verbosity verbosity)
{
  Scoped_Timer timer(timers, "Q");

  initialize_schur_off_diagonal(env, sdp, block_info, schur_complement,
                                schur_off_diagonal, schur_complement_cholesky,
                                timers, block_timings_ms, verbosity);
  syrk_Q(env, schur_off_diagonal, bigint_syrk_context, Q, timers,
         block_timings_ms, verbosity);
}