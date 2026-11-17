#include "catch2/catch_amalgamated.hpp"

#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/copy_matrix.hxx"
#include "unit_tests/util/util.hxx"

#include <El.hpp>

#include <memory>
#include <vector>

// cholesky_solve and lower_triangular_inverse_congruence
// (Abstract_Block_Diagonal_Matrix.hxx) call El::Trsm with TRSM_LARGE instead
// of Elemental's default, which is TRSM_MEDIUM when the right-hand side is no
// wider than 5 * (grid size). The result must be bitwise the same, on any grid
// shape and on both sides of that threshold.

namespace
{
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

  bool equal(const El::DistMatrix<El::BigFloat> &A,
             const El::DistMatrix<El::BigFloat> &B)
  {
    int eq = A.LocalHeight() == B.LocalHeight()
             && A.LocalWidth() == B.LocalWidth();
    for(El::Int j = 0; eq && j < A.LocalWidth(); ++j)
      for(El::Int i = 0; i < A.LocalHeight(); ++i)
        if(A.GetLocalCRef(i, j) != B.GetLocalCRef(i, j))
          eq = 0;
    return El::mpi::AllReduce(eq, El::mpi::MIN, El::mpi::COMM_WORLD) == 1;
  }
}

TEST_CASE("TRSM_LARGE in cholesky_solve and the step-length congruence: "
          "bitwise as El's default",
          "[trsm_large]")
{
  const int size = El::mpi::Size(El::mpi::COMM_WORLD);
  // size x 1, 1 x size and the default grid
  const int height = GENERATE(0, 1, -1);
  // below, at and above 5 * size (Medium / Large in El's default), and
  // 300: more than two of El's 128-row panels
  const int n
    = GENERATE_COPY(3, 2 * size + 1, 5 * size, 5 * size + 1, 37, 300);
  // El's block size: the default (128) and 8, so that small n also spans
  // several panels
  const int block_size = GENERATE(0, 8);
  DYNAMIC_SECTION("n=" << n << " height=" << height
                       << " block_size=" << block_size)
  {
    // restored even if a REQUIRE below fails
    struct Blocksize_Guard
    {
      int old = El::Blocksize();
      ~Blocksize_Guard() { El::SetBlocksize(old); }
    } guard;
    if(block_size > 0)
      El::SetBlocksize(block_size);
    const int h = height == -1 ? size : height;
    El::Grid grid(El::mpi::COMM_WORLD, h > 0 ? h : El::Grid::DefaultHeight(size));
    const auto A = same_matrix(n, true), R = same_matrix(n, false);
    Block_Diagonal_Matrix L({size_t(n)}, {0}, grid), X(L), Y(L), C(L);
    copy_matrix(A, L.blocks[0]);
    El::Cholesky(El::LOWER, L.blocks[0]);
    El::MakeTrapezoidal(El::LOWER, L.blocks[0]);
    copy_matrix(R, X.blocks[0]);
    Y = X;
    C = X;

    cholesky_solve(L, X);
    El::cholesky::SolveAfter(El::LOWER, El::NORMAL, L.blocks[0], Y.blocks[0]);
    REQUIRE(equal(X.blocks[0], Y.blocks[0]));

    auto D = C;
    lower_triangular_inverse_congruence(L, C);
    El::Trsm(El::RIGHT, El::LOWER, El::TRANSPOSE, El::NON_UNIT,
             El::BigFloat(1), L.blocks[0], D.blocks[0]);
    El::Trsm(El::LEFT, El::LOWER, El::NORMAL, El::NON_UNIT, El::BigFloat(1),
             L.blocks[0], D.blocks[0]);
    REQUIRE(equal(C.blocks[0], D.blocks[0]));
  }
}
