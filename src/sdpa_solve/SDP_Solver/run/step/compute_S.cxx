#include "sdp_solve/Block_Matrix/Block_Matrix.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpa_solve/SDP.hxx"
#include "sdpa_solve/memory_estimates.hxx"
#include "sdpb_util/Timers/Timers.hxx"

// Q = P^T P
// See sdp_solve/SDP_Solver/run/step/initialize_schur_complement_solver/compute_Q.cxx
// TODO rename
void syrk_Q(const Environment &env, Block_Matrix &P,
            BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
            El::DistMatrix<El::BigFloat> &Q, Timers &timers,
            El::Matrix<int32_t> &block_timings_ms, const Verbosity verbosity);

namespace Sdpb::Sdpa
{
  // compute P = (vec(G_1) vec(G_2) vec(G_3) ... )
  // G_p = L_X_Inv F_p L_Y
  // P is split vertically into blocks
  Block_Matrix
  initialize_P(const Environment &env, const SDP &sdp,
               const Block_Info &block_info, const El::Grid &grid,
               const Block_Diagonal_Matrix &X_cholesky,
               const Block_Diagonal_Matrix &Y_cholesky, Timers &timers,
               El::Matrix<int32_t> &block_timings_ms,
               const Verbosity verbosity)
  {
    Scoped_Timer timer(timers, "initialize_P");
    const size_t primal_dimension = block_info.primal_dimension;

    // Set matrix dimensions for output Block_Matrix P
    // Each block is a (dim*dim) x primal_dimension matrix
    std::vector<size_t> P_block_heights;
    P_block_heights.reserve(block_info.num_blocks());
    for(const auto &dim : block_info.block_dimensions)
      {
        P_block_heights.emplace_back(dim * dim);
      }
    const size_t P_width = primal_dimension;
    Block_Matrix result(P_block_heights, P_width, block_info.block_indices,
                        grid);
    if(verbosity >= Verbosity::trace)
      {
        print_allocation_message_per_node(env, "P",
                                          get_allocated_bytes(result));
      }

    for(size_t block = 0; block < block_info.num_blocks_local(); ++block)
      {
        const auto global_block_index = block_info.block_indices.at(block);
        auto block_index_string = std::to_string(global_block_index);
        const size_t block_dim
          = block_info.block_dimensions.at(global_block_index);

        Scoped_Timer solve_timer(timers, "solve_" + block_index_string);

        // Compute G_p = L_X_Inv F_p L_Y
        // start with F_p

        // Take blocks with current block index from F_1,F_2,..F_m
        // and put them on top of each other.
        El::DistMatrix<El::BigFloat> FY_block_vertical(
          block_dim * primal_dimension, block_dim, grid);
        for(size_t p = 0; p < primal_dimension; ++p)
          {
            const auto &source = sdp.sdp_blocks_F.at(p).blocks.at(block);
            const El::Range<int> I(p * block_dim, (p + 1) * block_dim);
            const El::Range<int> J(0, block_dim);
            auto dest = El::View(FY_block_vertical, I, J);
            dest = source;
            // sanity check
            ASSERT(dest.Viewing());
          }

        // compute F_p L_Y
        {
          Scoped_Timer trmm_timer(timers, "trmm");
          El::Trmm(El::LeftOrRightNS::RIGHT, El::UpperOrLowerNS::LOWER,
                   El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
                   El::BigFloat(1), Y_cholesky.blocks.at(block),
                   FY_block_vertical);
        }

        // compute L_X_inv F_p L_Y

        // Arrange blocks horizontally instead of vertically.
        // NB: ve do not transpose each block itself!
        Scoped_Timer transpose_G_timer(timers, "transpose_G");
        El::DistMatrix<El::BigFloat> G_block_horizontal(
          block_dim, block_dim * primal_dimension, grid);
        for(size_t p = 0; p < primal_dimension; ++p)
          {
            const El::Range<int> I(p * block_dim, (p + 1) * block_dim);
            const El::Range<int> J(0, block_dim);

            const auto source = FY_block_vertical(I, J);
            // NB: I <-> J
            auto dest = El::View(G_block_horizontal, J, I);
            dest = source;
            //sanity check
            ASSERT(dest.Viewing());
          }
        transpose_G_timer.stop();

        {
          Scoped_Timer trsm_timer(timers, "trsm");

          El::Trsm(El::LeftOrRightNS::LEFT, El::UpperOrLowerNS::LOWER,
                   El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
                   El::BigFloat(1), X_cholesky.blocks.at(block),
                   G_block_horizontal);
        }

        Scoped_Timer reshape_timer(timers, "reshape");
        for(size_t p = 0; p < sdp.primal_dimension(); ++p)
          {
            const El::Range<int> in_I(0, block_dim);
            const El::Range<int> in_J(p * block_dim, (p + 1) * block_dim);
            const auto G_p_block = El::View(G_block_horizontal, in_I, in_J);

            auto &P_block = result.blocks.at(block);
            // Write all elements of G_p into the p-th column of the resulting P matrix
            const int vec_height = P_block.Height();
            constexpr int vec_width = 1;
            const auto out_I = El::Range<int>(0, vec_height);
            const auto out_J = El::Range<int>(p, p + vec_width);

            auto vec_G_p_view = El::View(P_block, out_I, out_J);
            // TODO how does the resulting DistMatrix grid change after reshape?
            // Will the elements be evenly distributed among ranks?
            El::Reshape(vec_height, vec_width, G_p_block, vec_G_p_view);
          }
        block_timings_ms(global_block_index, 0)
          += solve_timer.elapsed_milliseconds();
      }

    return result;
  }

  void compute_S(const Environment &env, const SDP &sdp,
                 const Block_Info &block_info,
                 const Block_Diagonal_Matrix &X_cholesky,
                 const Block_Diagonal_Matrix &Y_cholesky, const El::Grid &grid,
                 BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
                 El::DistMatrix<El::BigFloat> &S, Timers &timers,
                 El::Matrix<int32_t> &block_timings_ms, Verbosity verbosity)
  {
    Scoped_Timer timer(timers, "S");
    Block_Matrix P
      = initialize_P(env, sdp, block_info, grid, X_cholesky, Y_cholesky,
                     timers, block_timings_ms, verbosity);
    syrk_Q(env, P, bigint_syrk_context, S, timers, block_timings_ms,
           verbosity);
  }
}
