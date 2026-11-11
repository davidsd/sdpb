#include "sdp_solve/SDP.hxx"
#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/Timers/Timers.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/memory_estimates.hxx"
#include "sdp_solve/SDP_Solver/run/multimodular/multimodular_blas.hxx"

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

void compute_schur_complement(
  const Block_Info &block_info,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  Block_Diagonal_Matrix &schur_complement, Timers &timers);

void compute_Q(const Environment &env, const SDP &sdp,
               const Block_Info &block_info,
               const Block_Diagonal_Matrix &schur_complement,
               Block_Matrix &schur_off_diagonal,
               Block_Diagonal_Matrix &schur_complement_cholesky,
               BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
               El::DistMatrix<El::BigFloat> &Q, Timers &timers,
               El::Matrix<int32_t> &block_timings_ms, bool multimodular_trsm,
               Verbosity verbosity);

void initialize_schur_complement_solver(
  const Environment &env, const Block_Info &block_info, const SDP &sdp,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const El::Grid &group_grid, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal,
  BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
  El::DistMatrix<El::BigFloat> &Q, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms, const bool multimodular_trsm,
  const int64_t multimodular_cholesky_Q, const Verbosity verbosity)
{
  Scoped_Timer initialize_timer(timers, "initializeSchurComplementSolver");
  // The Schur complement matrix S: a Block_Diagonal_Matrix with one
  // block for each 0 <= j < J.  SchurComplement.blocks[j] has dimension
  // (d_j+1)*m_j*(m_j+1)/2
  //
  Block_Diagonal_Matrix schur_complement(block_info.schur_block_sizes(),
                                         block_info.block_indices, group_grid);
  if(verbosity >= Verbosity::trace)
    {
      print_allocation_message_per_node(env, "schur_complement",
                                        get_allocated_bytes(schur_complement));
    }
  compute_schur_complement(block_info, A_X_inv, A_Y, schur_complement, timers);

  compute_Q(env, sdp, block_info, schur_complement, schur_off_diagonal,
            schur_complement_cholesky, bigint_syrk_context, Q, timers,
            block_timings_ms, multimodular_trsm, verbosity);

  Scoped_Timer Cholesky_timer(timers, "Cholesky_Q");
  try
    {
      if(multimodular_cholesky_Q > 0)
        {
          // El's blocked Cholesky with the trailing updates and the
          // diagonal blocks by multimodular BLAS, see
          // multimodular_cholesky_dist.cxx
          Multimodular_Blas_Config config;
          apply_multimodular_test_overrides(config);
          Multimodular_Cholesky_Dist_Timings t;
          // No per-panel timers (they would add ~5 entries per panel per
          // iteration to the timer list); Cholesky_Q times the whole.
          multimodular_cholesky_upper_dist(Q, int(multimodular_cholesky_Q),
                                           true, config, t, nullptr);
        }
      else
        Cholesky(El::UpperOrLowerNS::UPPER, Q);
    }
  catch(std::exception &e)
    {
      RUNTIME_ERROR("Error when computing Cholesky(Q): ", e.what());
    }
}
