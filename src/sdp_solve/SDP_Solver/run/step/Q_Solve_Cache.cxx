#include "Q_Solve_Cache.hxx"

#include "sdp_solve/SDP_Solver/run/multimodular/Node_Diag_Exchange.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <memory>
#include <optional>

using namespace El;
using F = BigFloat;

namespace
{
  // El's trsv::UT (src/blas_like/level2/Trsv/UT.hpp), x.Width() == 1
  // branch, with U11_STAR_STAR = U11 replaced by the kept copy.
  void trsv_ut(const Orientation orientation, const DistMatrix<F> &U,
               DistMatrix<F> &x, const std::vector<Matrix<F>> &blocks)
  {
    const Int m = U.Height();
    const Int bsize = Blocksize();
    const Grid &g = U.Grid();

    DistMatrix<F> U12(g), x1(g);
    DistMatrix<F, STAR, STAR> x1_STAR_STAR(g);

    DistMatrix<F, MC, STAR> x1_MC_STAR(g);
    DistMatrix<F, MR, MC> z1_MR_MC(g);
    DistMatrix<F, MR, STAR> z_MR_STAR(g);
    DistMatrix<F, MR, STAR> z1_MR_STAR(g), z2_MR_STAR(g);

    z_MR_STAR.AlignWith(U);
    z_MR_STAR.Resize(m, 1);
    Zero(z_MR_STAR);

    for(Int k = 0; k < m; k += bsize)
      {
        const Int nb = Min(bsize, m - k);

        LockedView(U12, U, IR(k, k + nb), IR(k + nb, m));

        View(x1, x, IR(k, k + nb), ALL);

        View(z1_MR_STAR, z_MR_STAR, IR(k, k + nb), ALL);
        View(z2_MR_STAR, z_MR_STAR, IR(k + nb, m), ALL);

        if(k != 0)
          {
            Contract(z1_MR_STAR, z1_MR_MC);
            x1 += z1_MR_MC;
          }

        x1_STAR_STAR = x1;
        Trsv(UPPER, orientation, NON_UNIT, blocks.at(k / bsize),
             x1_STAR_STAR.Matrix());
        x1 = x1_STAR_STAR;

        x1_MC_STAR.AlignWith(U12);
        x1_MC_STAR = x1_STAR_STAR;
        LocalGemv(orientation, F(-1), U12, x1_MC_STAR, F(1), z2_MR_STAR);
      }
  }

  // El's trsv::UN (src/blas_like/level2/Trsv/UN.hpp), x.Width() == 1
  // branch, with U11_STAR_STAR = U11 replaced by the kept copy.
  void trsv_un(const DistMatrix<F> &U, DistMatrix<F> &x,
               const std::vector<Matrix<F>> &blocks)
  {
    const Int m = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset(m, bsize);
    const Grid &g = U.Grid();

    DistMatrix<F> U01(g), x1(g);
    DistMatrix<F, STAR, STAR> x1_STAR_STAR(g);

    DistMatrix<F, MR, STAR> x1_MR_STAR(g);
    DistMatrix<F, MC, STAR> z_MC_STAR(g);
    DistMatrix<F, MC, STAR> z0_MC_STAR(g), z1_MC_STAR(g);

    z_MC_STAR.AlignWith(U);
    z_MC_STAR.Resize(m, 1);
    Zero(z_MC_STAR);

    for(Int k = kLast; k >= 0; k -= bsize)
      {
        const Int nb = Min(bsize, m - k);

        LockedView(U01, U, IR(0, k), IR(k, k + nb));

        View(x1, x, IR(k, k + nb), ALL);

        View(z0_MC_STAR, z_MC_STAR, IR(0, k), ALL);
        View(z1_MC_STAR, z_MC_STAR, IR(k, k + nb), ALL);

        if(k + nb != m)
          AxpyContract(F(1), z1_MC_STAR, x1);

        x1_STAR_STAR = x1;
        Trsv(UPPER, NORMAL, NON_UNIT, blocks.at(k / bsize),
             x1_STAR_STAR.Matrix());
        x1 = x1_STAR_STAR;

        x1_MR_STAR.AlignWith(U01);
        x1_MR_STAR = x1_STAR_STAR;
        LocalGemv(NORMAL, F(-1), U01, x1_MR_STAR, F(1), z0_MC_STAR);
      }
  }
}

Q_Solve_Cache::Q_Solve_Cache(const bool node_shared_copies)
    : node_shared(node_shared_copies)
{}

void Q_Solve_Cache::build(const DistMatrix<F> &U, Timers *timers)
{
  std::optional<Scoped_Timer> timer;
  if(timers != nullptr)
    timer.emplace(*timers, "Q_solve_cache_build");
  ASSERT(U.ColDist() == MC && U.RowDist() == MR,
         "Q_Solve_Cache: expected an [MC,MR] matrix");
  const Int m = U.Height();
  const Int bsize = Blocksize();
  blocks.clear();
  blocks.reserve((m + bsize - 1) / bsize);
  std::unique_ptr<Node_Diag_Exchange> exchange;
  if(node_shared && m > 0)
    exchange = std::make_unique<Node_Diag_Exchange>(
      U.Grid(), size_t(Min(bsize, m)) * size_t(Min(bsize, m)));
  DistMatrix<F> U11(U.Grid());
  DistMatrix<F, STAR, STAR> U11_STAR_STAR(U.Grid());
  for(Int k = 0; k < m; k += bsize)
    {
      const Int nb = Min(bsize, m - k);
      LockedView(U11, U, IR(k, k + nb), IR(k, k + nb));
      if(exchange)
        exchange->all_gather(U11, U11_STAR_STAR);
      else
        U11_STAR_STAR = U11;
      blocks.emplace_back(U11_STAR_STAR.LockedMatrix());
    }
  U_buffer = U.LockedBuffer();
  U_height = m;
}

void Q_Solve_Cache::solve(const DistMatrix<F> &U, DistMatrix<F> &x,
                          Timers *timers)
{
  ASSERT_EQUAL(x.Width(), 1);
  if(U_height != U.Height() || U_buffer != U.LockedBuffer())
    build(U, timers);
  // El::cholesky::SolveAfter(UPPER, NORMAL, U, x): Trsm(LEFT, UPPER,
  // ADJOINT, NON_UNIT, 1, U, x), then Trsm(LEFT, UPPER, NORMAL, ...); each
  // Trsm scales x by alpha = 1 and calls Trsv for one right-hand side.
  x *= F(1);
  trsv_ut(ADJOINT, U, x, blocks);
  x *= F(1);
  trsv_un(U, x, blocks);
}
