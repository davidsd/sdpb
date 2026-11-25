#include "catch2/catch_amalgamated.hpp"

#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/copy_matrix.hxx"
#include "test_util/diff.hxx"
#include "unit_tests/util/util.hxx"

#include <El.hpp>
#include <boost/noncopyable.hpp>

#include <vector>

using Test_Util::REQUIRE_Equal::diff;

// cholesky_solve and lower_triangular_inverse_congruence
// (Abstract_Block_Diagonal_Matrix.hxx) call El::Trsm with TRSM_LARGE instead
// of Elemental's default, which is TRSM_MEDIUM when the right-hand side is no
// wider than 5 * (grid size). The result must be bitwise the same as with
// TRSM_MEDIUM, on any grid shape. The reference asks for TRSM_MEDIUM
// explicitly, so that every case compares the two algorithms, whatever the
// number of ranks.

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

  // Bitwise comparison of two DistMatrices. Plain DIFF_PREC would fail only
  // on the ranks holding a mismatch, and the other ranks would then hang in
  // the next section's collectives. Here every rank fails together, and
  // DIFF_PREC reports the first mismatch on the ranks that have one.
  void diff_bitwise(const El::DistMatrix<El::BigFloat> &a,
                    const El::DistMatrix<El::BigFloat> &b)
  {
    // The local entries of a and b are the same global entries only if the
    // sizes, grids and alignments match. These are the same on every rank,
    // so a mismatch fails every rank.
    REQUIRE(a.Height() == b.Height());
    REQUIRE(a.Width() == b.Width());
    REQUIRE(a.Grid() == b.Grid());
    REQUIRE(a.ColAlign() == b.ColAlign());
    REQUIRE(a.RowAlign() == b.RowAlign());
    int local_equal = 1;
    for(El::Int jLoc = 0; local_equal && jLoc < a.LocalWidth(); ++jLoc)
      for(El::Int iLoc = 0; iLoc < a.LocalHeight(); ++iLoc)
        if(a.GetLocalCRef(iLoc, jLoc) != b.GetLocalCRef(iLoc, jLoc))
          local_equal = 0;
    const bool all_equal
      = El::mpi::AllReduce(local_equal, El::mpi::MIN, a.DistComm()) == 1;
    if(!all_equal)
      DIFF_PREC(a, b, -1);
    REQUIRE(all_equal);
  }

  // Restores El's block size, even if a REQUIRE fails
  struct Blocksize_Guard : boost::noncopyable
  {
    const int old = El::Blocksize();
    ~Blocksize_Guard() { El::SetBlocksize(old); }
  };
}

TEST_CASE("TRSM_LARGE in cholesky_solve and the step-length congruence: "
          "bitwise as TRSM_MEDIUM",
          "[trsm_large]")
{
  const int size = El::mpi::Size(El::mpi::COMM_WORLD);
  // Grid height: El's default, 1 (1 x size), or size (size x 1)
  const int height = GENERATE_COPY(0, 1, size);
  const int n = GENERATE(3, 13, 37);
  // El's block size: the default (128) and 8, so that n also spans several
  // panels
  const int block_size = GENERATE(0, 8);
  DYNAMIC_SECTION("n=" << n << " height=" << height
                       << " block_size=" << block_size)
  {
    Blocksize_Guard guard;
    if(block_size > 0)
      El::SetBlocksize(block_size);
    El::Grid grid(El::mpi::COMM_WORLD,
                  height > 0 ? height : El::Grid::DefaultHeight(size));
    const auto A = same_matrix(n, true), R = same_matrix(n, false);
    Block_Diagonal_Matrix L({size_t(n)}, {0}, grid), X(L), Y(L), C(L);
    copy_matrix(A, L.blocks[0]);
    El::Cholesky(El::LOWER, L.blocks[0]);
    El::MakeTrapezoidal(El::LOWER, L.blocks[0]);
    copy_matrix(R, X.blocks[0]);
    Y = X;
    C = X;

    // El::cholesky::SolveAfter(LOWER, NORMAL), with TRSM_MEDIUM
    cholesky_solve(L, X);
    El::Trsm(El::LEFT, El::LOWER, El::NORMAL, El::NON_UNIT, El::BigFloat(1),
             L.blocks[0], Y.blocks[0], false, El::TRSM_MEDIUM);
    El::Trsm(El::LEFT, El::LOWER, El::ADJOINT, El::NON_UNIT, El::BigFloat(1),
             L.blocks[0], Y.blocks[0], false, El::TRSM_MEDIUM);
    diff_bitwise(X.blocks[0], Y.blocks[0]);

    auto D = C;
    lower_triangular_inverse_congruence(L, C);
    El::Trsm(El::RIGHT, El::LOWER, El::TRANSPOSE, El::NON_UNIT,
             El::BigFloat(1), L.blocks[0], D.blocks[0]);
    El::Trsm(El::LEFT, El::LOWER, El::NORMAL, El::NON_UNIT, El::BigFloat(1),
             L.blocks[0], D.blocks[0], false, El::TRSM_MEDIUM);
    diff_bitwise(C.blocks[0], D.blocks[0]);
  }
}
