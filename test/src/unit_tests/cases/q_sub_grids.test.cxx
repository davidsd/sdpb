#include "catch2/catch_amalgamated.hpp"

#include "sdp_solve/SDP_Solver/run/multimodular/multimodular_blas.hxx"
#include "sdp_solve/SDP_Solver/run/step/Q_Grids.hxx"
#include "sdpb_util/copy_matrix.hxx"
#include "unit_tests/util/util.hxx"

#include <El.hpp>

#include <cmath>
#include <memory>
#include <vector>

// SDPB's Q on a sub-grid of the ranks (--qSolveRanks): copies between
// grids, and factorizations that don't depend on the grid.
namespace
{
  // Grid on world ranks first, ..., first + r - 1 (mod size), viewed by all
  // ranks, h x (r / h) as in Q_Grids, or height x (r / height)
  std::unique_ptr<El::Grid>
  sub_grid(const int first, const int r, const int height_in = 0)
  {
    const int size = El::mpi::Size(El::mpi::COMM_WORLD);
    std::vector<int> ranks(r);
    for(int i = 0; i < r; ++i)
      ranks[i] = (first + i) % size;
    int height = int(std::sqrt(double(r)));
    while(r % height != 0)
      --height;
    if(height_in > 0)
      height = height_in;
    El::mpi::Group world, owners;
    El::mpi::CommGroup(El::mpi::COMM_WORLD, world);
    El::mpi::Incl(world, r, ranks.data(), owners);
    auto grid
      = std::make_unique<El::Grid>(El::mpi::COMM_WORLD, owners, height);
    El::mpi::Free(owners);
    El::mpi::Free(world);
    return grid;
  }

  // The same matrix on every rank: random (spd: G G^T + n I) on rank 0,
  // broadcast
  El::Matrix<El::BigFloat> same_matrix(const int n, const bool spd)
  {
    El::Matrix<El::BigFloat> M(n, n);
    if(El::mpi::Rank() == 0)
      {
        const auto G = Test_Util::random_matrix(n, n);
        if(spd)
          {
            El::Gemm(El::NORMAL, El::TRANSPOSE, El::BigFloat(1), G, G,
                     El::BigFloat(0), M);
            for(int i = 0; i < n; ++i)
              M(i, i) += El::BigFloat(n);
          }
        else
          M = G;
      }
    El::mpi::Broadcast(M.Buffer(), n * n, 0, El::mpi::COMM_WORLD);
    return M;
  }

  // Every local entry of B: the upper triangle of M, zero below
  bool is_upper_of(const El::DistMatrix<El::BigFloat> &B,
                   const El::Matrix<El::BigFloat> &M)
  {
    for(El::Int jLoc = 0; jLoc < B.LocalWidth(); ++jLoc)
      for(El::Int iLoc = 0; iLoc < B.LocalHeight(); ++iLoc)
        {
          const El::Int i = B.GlobalRow(iLoc), j = B.GlobalCol(jLoc);
          const El::BigFloat want = i <= j ? M.CRef(i, j) : El::BigFloat(0);
          if(B.GetLocalCRef(iLoc, jLoc) != want)
            return false;
        }
    return true;
  }

  // Grids of a test: sub-grids of 1, size/2 and size-1 ranks starting at
  // rank 0 or 1, and all ranks as size x 1 and 1 x size
  std::vector<std::unique_ptr<El::Grid>> test_grids()
  {
    const int size = El::mpi::Size(El::mpi::COMM_WORLD);
    std::vector<std::unique_ptr<El::Grid>> grids;
    for(const int r : {1, size / 2, size - 1})
      for(const int first : {0, 1})
        if(r >= 1 && first < size)
          grids.push_back(sub_grid(first, r));
    if(size > 1)
      {
        grids.push_back(sub_grid(0, size, size));
        grids.push_back(sub_grid(0, size, 1));
      }
    return grids;
  }
}

TEST_CASE("copy_upper_between_grids", "[q_sub_grids]")
{
  const int n = GENERATE(1, 5, 37, 130);
  DYNAMIC_SECTION("n=" << n)
  {
    const auto M = same_matrix(n, false);
    const auto grids = test_grids();
    std::vector<const El::Grid *> all{&El::Grid::Default()};
    for(const auto &g : grids)
      all.push_back(g.get());
    for(size_t a = 0; a < all.size(); ++a)
      for(size_t b = 0; b < all.size(); ++b)
        {
          CAPTURE(a, b, all[a]->Size(), all[b]->Size());
          El::DistMatrix<El::BigFloat> A(*all[a]), B(*all[b]);
          A.Resize(n, n);
          copy_matrix(M, A);
          copy_upper_between_grids(A, B);
          REQUIRE(B.Height() == n);
          REQUIRE(B.Width() == n);
          REQUIRE(is_upper_of(B, M));
        }
  }
}

// Q's Cholesky on a sub-grid is a pure layout change: El::Cholesky(UPPER)
// (variant 3: diagonal blocks on [STAR,STAR], panel solves per column,
// trailing updates of length nb per entry) and the multimodular Cholesky_Q
// (row/column scales over the panel width) give bitwise the same factor on
// any grid.
TEST_CASE("Cholesky(Q) on a sub-grid: same factor", "[q_sub_grids]")
{
  const int n = GENERATE(37, 200);
  // 0: El::Cholesky; else multimodular, panel width, diagonal blocks by El
  // (+0) or multimodular (+1)
  const int variant = GENERATE(0, 32, 33, 64);
  DYNAMIC_SECTION("n=" << n << " variant=" << variant)
  {
    const auto M = same_matrix(n, true);
    auto factor = [&](El::DistMatrix<El::BigFloat> &U) {
      if(!U.Participating())
        return;
      if(variant == 0)
        El::Cholesky(El::UPPER, U);
      else
        {
          Multimodular_Blas_Config config;
          Multimodular_Cholesky_Dist_Timings timings;
          multimodular_cholesky_upper_dist(U, variant & ~1, variant & 1,
                                           config, timings);
        }
    };
    El::DistMatrix<El::BigFloat> U_all(El::Grid::Default());
    U_all.Resize(n, n);
    copy_matrix(M, U_all);
    factor(U_all);
    El::Matrix<El::BigFloat> U_ref(n, n);
    El::Zero(U_ref);
    {
      const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> U_star(U_all);
      for(int j = 0; j < n; ++j)
        for(int i = 0; i <= j; ++i)
          U_ref(i, j) = U_star.LockedMatrix().CRef(i, j);
    }
    for(const auto &g : test_grids())
      {
        CAPTURE(g->Size());
        El::DistMatrix<El::BigFloat> U(*g);
        U.Resize(n, n);
        copy_matrix(M, U);
        factor(U);
        // back to all ranks, upper triangle only
        El::DistMatrix<El::BigFloat> U_back(El::Grid::Default());
        copy_upper_between_grids(U, U_back);
        REQUIRE(is_upper_of(U_back, U_ref));
      }
  }
}
