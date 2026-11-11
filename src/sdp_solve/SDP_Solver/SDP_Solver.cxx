#include "../SDP_Solver.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Boost_Float.hxx"

#include <tuple>

namespace
{
  // --initialSeed=rows[:a]. In each Gram
  // sector (block j, parity) of a block with n > 1 samples,
  //
  //   X = Omega_p 1_m (x) M,  Y = Omega_d 1_m (x) M^-1,
  //   M = 1 + Q diag(g^a - 1) Q^T,
  //
  // Q the sector's bilinear basis (d x n), m the block's matrix dimension,
  // and g_k = 1 / max(|c_p|, |B_pn|) over the block's rows p = (r,s,k) at
  // sample k. pmp2sdp's Q has orthonormal rows, so M = Q diag(g^a) Q^T to
  // that accuracy, and a = 0 gives M = 1 exactly: the default start. A
  // one-sample block keeps M = 1 (pmp2sdp fixes its basis to 1).
  // XY = Omega_p Omega_d 1, so mu starts where the default start's does;
  // x = y = 0.
  void apply_rows_seed(const std::string &seed, const Solver_Parameters &p,
                       const Block_Info &block_info, const El::Grid &grid,
                       const SDP &sdp, Paired_Block_Diagonal_Matrix &X,
                       Paired_Block_Diagonal_Matrix &Y,
                       const Verbosity &verbosity)
  {
    const Boost_Float a(seed.size() > 5 ? seed.substr(5) : "0.5");
    const auto &indices(block_info.block_indices);

    // Largest |c_p|, |B_pn| at each sample of each local block; one
    // reduction over the group of ranks that holds these blocks.
    std::vector<size_t> offset(indices.size() + 1, 0);
    for(size_t b = 0; b < indices.size(); ++b)
      offset[b + 1] = offset[b] + block_info.num_points.at(indices[b]);
    std::vector<El::BigFloat> largest(offset.back(), El::BigFloat(0));
    for(size_t b = 0; b < indices.size(); ++b)
      {
        const int64_t n(block_info.num_points.at(indices[b]));
        for(auto *rows : {&sdp.primal_objective_c.blocks.at(b),
                          &sdp.free_var_matrix.blocks.at(b)})
          for(int64_t r = 0; r < rows->LocalHeight(); ++r)
            for(int64_t c = 0; c < rows->LocalWidth(); ++c)
              {
                auto &l(largest[offset[b] + rows->GlobalRow(r) % n]);
                l = std::max(l, El::BigFloat(El::Abs(rows->GetLocal(r, c))));
              }
      }
    if(!largest.empty())
      El::mpi::AllReduce(largest.data(), int(largest.size()), El::mpi::MAX,
                         grid.Comm());
    auto g = [&](const size_t b, const size_t k) {
      const El::BigFloat &l(largest[offset[b] + k]);
      return l > El::BigFloat(0) ? El::BigFloat(1) / l : El::BigFloat(1);
    };

    El::BigFloat g_min(1), g_max(1);
    for(size_t b = 0; b < indices.size(); ++b)
      {
        const int64_t n(block_info.num_points.at(indices[b])),
          m(block_info.dimensions.at(indices[b]));
        if(n < 2)
          continue;
        for(int64_t k = 0; k < n; ++k)
          {
            g_min = std::min(g_min, g(b, k));
            g_max = std::max(g_max, g(b, k));
          }
        for(size_t parity = 0; parity < 2; ++parity)
          {
            const auto &Q(sdp.bilinear_bases.at(2 * b + parity));
            const int64_t d(Q.Height());
            if(d == 0)
              continue;
            ASSERT_EQUAL(Q.Width(), n);
            // M is positive definite (= Q diag(g^a) Q^T) only if Q has
            // orthonormal rows, as pmp2sdp's default bases do; a PMP may give
            // its own bilinearBasis.
            {
              El::DistMatrix<El::BigFloat> QQT(d, d, grid);
              El::Gemm(El::Orientation::NORMAL, El::Orientation::TRANSPOSE,
                       El::BigFloat(1), Q, Q, El::BigFloat(0), QQT);
              El::ShiftDiagonal(QQT, El::BigFloat(-1));
              // pmp2sdp's bases are orthonormal to working precision
              El::BigFloat tolerance(1);
              mpf_div_2exp(tolerance.gmp_float.get_mpf_t(),
                           tolerance.gmp_float.get_mpf_t(),
                           El::gmp::Precision() / 2);
              if(El::MaxAbs(QQT) > tolerance)
                RUNTIME_ERROR(
                  "--initialSeed=", seed, " needs bilinear bases with "
                  "orthonormal rows (pmp2sdp's default ones); block ",
                  indices[b], " has |Q Q^T - 1| = ", El::MaxAbs(QQT), ".");
            }
            El::DistMatrix<El::BigFloat> w(n, 1, grid);
            for(int64_t r = 0; r < w.LocalHeight(); ++r)
              for(int64_t c = 0; c < w.LocalWidth(); ++c)
                w.SetLocal(r, c,
                           to_BigFloat(pow(to_Boost_Float(g(b, w.GlobalRow(r))),
                                           a))
                             - El::BigFloat(1));
            El::DistMatrix<El::BigFloat> scaled(Q), M(d, d, grid);
            El::DiagonalScale(El::LeftOrRight::RIGHT, El::Orientation::NORMAL,
                              w, scaled);
            El::Gemm(El::Orientation::NORMAL, El::Orientation::TRANSPOSE,
                     El::BigFloat(1), Q, scaled, El::BigFloat(0), M);
            El::ShiftDiagonal(M, El::BigFloat(1));
            // M^-1 = L^-T L^-1 with M = L L^T.
            El::DistMatrix<El::BigFloat> L(M), M_inv(grid);
            El::Cholesky(El::UpperOrLowerNS::LOWER, L);
            El::Identity(M_inv, d, d);
            El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
                     El::OrientationNS::NORMAL, El::UnitOrNonUnit::NON_UNIT,
                     El::BigFloat(1), L, M_inv);
            El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
                     El::OrientationNS::TRANSPOSE, El::UnitOrNonUnit::NON_UNIT,
                     El::BigFloat(1), L, M_inv);
            // Each rank sets its own entries of X and Y from the lower
            // triangles, so both are exactly symmetric.
            const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> M_star(M),
              M_inv_star(M_inv);
            for(auto [S, Z, scale] :
                {std::make_tuple(&X.get_block(b, parity), &M_star,
                                 p.initial_matrix_scale_primal),
                 std::make_tuple(&Y.get_block(b, parity), &M_inv_star,
                                 p.initial_matrix_scale_dual)})
              {
                ASSERT_EQUAL(S->Height(), m * d);
                for(int64_t r = 0; r < S->LocalHeight(); ++r)
                  for(int64_t c = 0; c < S->LocalWidth(); ++c)
                    {
                      const int64_t R(S->GlobalRow(r)), C(S->GlobalCol(c));
                      S->SetLocal(r, c,
                                  R / d != C / d
                                    ? El::BigFloat(0)
                                    : scale
                                        * Z->GetLocal(std::max(R % d, C % d),
                                                      std::min(R % d, C % d)));
                    }
              }
          }
      }
    if(verbosity >= Verbosity::regular)
      {
        const El::BigFloat lo(
          El::mpi::AllReduce(g_min, El::mpi::MIN, El::mpi::COMM_WORLD)),
          hi(El::mpi::AllReduce(g_max, El::mpi::MAX, El::mpi::COMM_WORLD));
        if(El::mpi::Rank() == 0)
          El::Output("Initial seed ", seed,
                     ": X = Omega_p M, Y = Omega_d M^-1, "
                     "M = 1 + Q diag(g^a - 1) Q^T, a = ",
                     to_string(a), ", g from the SDP's rows, ", lo, " to ", hi,
                     " (1 included)");
      }
  }
}

SDP_Solver::SDP_Solver(const Solver_Parameters &parameters,
                       const Verbosity &verbosity,
                       const bool &require_initial_checkpoint,
                       const Block_Info &block_info, const El::Grid &grid,
                       const size_t &dual_objective_b_height, const SDP *sdp)
    : x(block_info.schur_block_sizes(), block_info.block_indices, grid),
      X(block_info.psd_matrix_block_sizes(), block_info.block_indices, grid),
      y(std::vector<size_t>(block_info.num_points.size(),
                            dual_objective_b_height),
        block_info.block_indices, grid),
      Y(X),
      primal_residues(X),
      dual_residues(block_info.schur_block_sizes(), block_info.block_indices,
                    grid),
      current_generation(0)
{
  if(!load_checkpoint(parameters.checkpoint_in, block_info, verbosity,
                      require_initial_checkpoint))
    {
      X.set_zero();
      Y.set_zero();
      for(auto &block : x.blocks)
        {
          El::Zero(block);
        }
      for(auto &block : y.blocks)
        {
          El::Zero(block);
        }

      // X = \Omega_p I
      X.add_diagonal(parameters.initial_matrix_scale_primal);
      // Y = \Omega_d I
      Y.add_diagonal(parameters.initial_matrix_scale_dual);

      // Opt-in seed: only on a cold start, never over a checkpoint.
      const std::string &seed(parameters.initial_seed);
      if(!seed.empty() && seed != "none")
        {
          if(sdp == nullptr)
            RUNTIME_ERROR("--initialSeed=", seed,
                          " needs the SDP, which this program doesn't pass");
          apply_rows_seed(seed, parameters, block_info, grid, *sdp, X, Y,
                          verbosity);
        }
    }
  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      const size_t num_blocks = block_info.num_points.size();
      size_t primal_dimension = 0;
      for(size_t i = 0; i < num_blocks; ++i)
        primal_dimension += block_info.get_schur_block_size(i);
      const size_t dual_dimension = dual_objective_b_height;
      El::Output("Initialize SDP solver"
                 "\n\tprimal dimension: ",
                 primal_dimension, "\n\tdual dimension: ", dual_dimension,
                 "\n\tSDP blocks: ", num_blocks);
    }
}
