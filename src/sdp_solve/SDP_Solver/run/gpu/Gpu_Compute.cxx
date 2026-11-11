#include "Gpu_Compute.hxx"
#include "gpu/rns/Limb_Sizes.hxx"
#include "gpu/solver/Gpu_Solver.hxx"
#include "sdp_solve/SDP_Solver.hxx"
#include "sdp_solve/SDP_Solver/run/step/centering_parameters.hxx"
#include "sdpb_util/assert.hxx"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <thread>
#include <vector>

using namespace sdpb_gpu;

namespace
{
  // Largest limb count the conversions below handle (their stack buffers)
  constexpr int max_limbs = 79;
}

int Gpu_Compute::limbs_for_precision(const size_t precision)
{
  // limbs for --precision, rounded up to a size with compiled kernels
  const int L = sdpb_gpu::rounded_L(int((precision + 2 + 31) / 32) + 1);
  ASSERT(L > 0, "--gpu: --precision ", precision,
         " is above the largest float size with device kernels (",
         sdpb_gpu::largest_compiled_L(),
         " limbs, gpu/rns/Limb_Sizes.hxx)");
  static_assert(sdpb_gpu::largest_compiled_L() <= max_limbs);
  return L;
}

namespace
{
  constexpr int64_t MPF_ZERO_EXP = -(int64_t(1) << 40);

  // v -> device float: exactly NB = 32 L - 2 bits, rounded half up on the magnitude
  inline void to_device(const El::BigFloat &v, int L, int32_t *m, int64_t plane,
                        int64_t idx, int64_t *e)
  {
    const __mpf_struct *f = v.gmp_float.get_mpf_t();
    const int size = f->_mp_size;
    if(size == 0)
      {
        for(int t = 0; t < L; ++t)
          m[t * plane + idx] = 0;
        e[idx] = MPF_ZERO_EXP;
        return;
      }
    const int n = size < 0 ? -size : size;
    const mp_limb_t *d = f->_mp_d;
    const int NB = 32 * L - 2;
    const int64_t bits = 64 * int64_t(n) - __builtin_clzll(d[n - 1]);
    int64_t sh = bits - NB;
    uint32_t r[max_limbs + 1];
    auto limb32 = [&](int64_t i) -> uint32_t {
      if(i < 0 || i >= 2 * int64_t(n))
        return 0;
      const mp_limb_t w = d[i / 2];
      return (i & 1) ? uint32_t(w >> 32) : uint32_t(w);
    };
    auto extract = [&](int64_t s) {
      for(int t = 0; t <= L; ++t)
        {
          const int64_t bit = 32 * int64_t(t) + s;
          const int64_t q = bit >= 0 ? bit / 32 : -((-bit + 31) / 32);
          const int rb = int(bit - q * 32);
          r[t] = uint32_t((((uint64_t)limb32(q + 1) << 32) | limb32(q)) >> rb);
        }
      if(s >= 1)
        {
          const int64_t gb = s - 1;
          uint64_t c = (limb32(gb / 32) >> (gb % 32)) & 1u;
          for(int t = 0; t <= L && c; ++t)
            {
              const uint64_t w = uint64_t(r[t]) + c;
              r[t] = uint32_t(w);
              c = w >> 32;
            }
        }
    };
    extract(sh);
    if(r[L] != 0 || (r[L - 1] >> (NB - 32 * (L - 1))) != 0)
      {
        ++sh;
        extract(sh);
      }
    if(size < 0)
      {
        uint64_t c = 1;
        for(int t = 0; t < L; ++t)
          {
            const uint64_t w = uint64_t(~r[t]) + c;
            r[t] = uint32_t(w);
            c = w >> 32;
          }
      }
    for(int t = 0; t < L; ++t)
      m[t * plane + idx] = int32_t(r[t]);
    // the integer of `bits` bits has exponent 64 (exp - n)
    e[idx] = 64 * (int64_t(f->_mp_exp) - n) + sh;
  }

  inline void one_device(int L, int32_t *m, int64_t plane, int64_t idx,
                         int64_t *e)
  {
    const int NB = 32 * L - 2;
    for(int t = 0; t < L; ++t)
      m[t * plane + idx] = t == L - 1 ? int32_t(1u << (NB - 1 - 32 * (L - 1))) : 0;
    e[idx] = -(NB - 1);
  }

  inline void zero_device(int L, int32_t *m, int64_t plane, int64_t idx,
                          int64_t *e)
  {
    for(int t = 0; t < L; ++t)
      m[t * plane + idx] = 0;
    e[idx] = MPF_ZERO_EXP;
  }

  // limbs (two's complement) * 2^e -> BigFloat, exact when it fits the precision
  El::BigFloat from_limbs(const uint32_t *limbs, int k, int64_t stride, int64_t e)
  {
    std::vector<uint32_t> a(k);
    for(int t = 0; t < k; ++t)
      a[t] = limbs[t * stride];
    const bool neg = int32_t(a[k - 1]) < 0;
    if(neg)
      {
        uint64_t c = 1;
        for(int t = 0; t < k; ++t)
          {
            const uint64_t w = uint64_t(~a[t]) + c;
            a[t] = uint32_t(w);
            c = w >> 32;
          }
      }
    mpz_t z;
    mpz_init(z);
    mpz_import(z, k, -1, 4, 0, 0, a.data());
    El::BigFloat r;
    mpf_set_z(r.gmp_float.get_mpf_t(), z);
    mpz_clear(z);
    if(e >= 0)
      mpf_mul_2exp(r.gmp_float.get_mpf_t(), r.gmp_float.get_mpf_t(), e);
    else
      mpf_div_2exp(r.gmp_float.get_mpf_t(), r.gmp_float.get_mpf_t(), -e);
    if(neg)
      mpf_neg(r.gmp_float.get_mpf_t(), r.gmp_float.get_mpf_t());
    return r;
  }

  El::BigFloat from_host(const Host_Float &h)
  {
    if(h.m.empty() || h.e == MPF_ZERO_EXP)
      return El::BigFloat(0);
    return from_limbs(h.m.data(), int(h.m.size()), 1, h.e);
  }

  Host_Float to_host(const El::BigFloat &v, int L)
  {
    Host_Float h;
    std::vector<int32_t> m(L);
    int64_t e;
    to_device(v, L, m.data(), 1, 0, &e);
    h.m.assign(reinterpret_cast<uint32_t *>(m.data()),
               reinterpret_cast<uint32_t *>(m.data()) + L);
    h.e = e;
    return h;
  }

  // fn(i) for 0 <= i < n on OMP_NUM_THREADS threads. Not OpenMP: Elemental's
  // headers change with _OPENMP, and they must read alike in all of SDPB.
  template <class Fn> void parallel_for(int64_t n, const Fn &fn)
  {
    const char *s = getenv("OMP_NUM_THREADS");
    int nt = s ? atoi(s) : int(std::thread::hardware_concurrency());
    nt = int(std::max<int64_t>(
      1, std::min<int64_t>(nt, std::max<int64_t>(1, n / 1024))));
    if(nt == 1)
      {
        for(int64_t i = 0; i < n; ++i)
          fn(i);
        return;
      }
    std::vector<std::thread> pool;
    for(int t = 0; t < nt; ++t)
      pool.emplace_back([&, t]() {
        for(int64_t i = n * t / nt; i < n * (t + 1) / nt; ++i)
          fn(i);
      });
    for(auto &th : pool)
      th.join();
  }

  long env_long(const char *name, long dflt)
  {
    const char *s = getenv(name);
    return s ? atol(s) : dflt;
  }
  double env_double(const char *name, double dflt)
  {
    const char *s = getenv(name);
    return s ? atof(s) : dflt;
  }
}

struct Gpu_Compute::Impl
{
  std::unique_ptr<Gpu_Solver> gs;
  int L, N;
  // local index of a global block (one rank: all blocks are local)
  std::vector<int> local_of;
  std::vector<int> dim_of;
  Residue_Result last;
  bool have_iterate_on_host = false;
  std::vector<size_t> block_dims;
  std::vector<int> he, ho, Kj;

  template <class Fn> void fill(int chunk, Gpu_Array a, const Fn &fn)
  {
    const int64_t tot = gs->entries(chunk, a);
    std::vector<int32_t> m(size_t(L) * tot);
    std::vector<int64_t> e(tot);
    parallel_for(tot, [&](int64_t q) { fn(q, m.data(), tot, e.data()); });
    gs->upload(chunk, a, m.data(), e.data());
  }

  // rows of c, B, x: segment s, point k -> s K + k
  static int block_row(int row, int K, int Kb)
  {
    const int s = row / K, k = row % K;
    return k < Kb ? s * Kb + k : -1;
  }

  void upload_b(const SDP &sdp)
  {
    std::vector<int32_t> m(size_t(L) * N);
    std::vector<int64_t> e(N);
    const auto &bm = sdp.dual_objective_b.LockedMatrix();
    for(int i = 0; i < N; ++i)
      to_device(bm.CRef(i, 0), L, m.data(), N, i, e.data());
    gs->upload_b(m.data(), e.data());
  }

  void upload_c_B(int ci, const SDP &sdp)
  {
    const Chunk_Info &I = gs->chunk(ci);
    const int K = I.K, Pc = I.P, Nv = N;
    fill(ci, Gpu_Array::c, [&](int64_t q, int32_t *m, int64_t tot, int64_t *e) {
      const int q_ = int(q / Pc), row = int(q % Pc);
      const int j = I.blocks[q_];
      const int r = block_row(row, K, Kj[j]);
      if(r >= 0)
        to_device(sdp.primal_objective_c.blocks[local_of[j]].LockedMatrix().CRef(r, 0),
                  L, m, tot, q, e);
      else
        zero_device(L, m, tot, q, e);
    });
    fill(ci, Gpu_Array::B, [&](int64_t q, int32_t *m, int64_t tot, int64_t *e) {
      const int64_t per = int64_t(Pc) * Nv;
      const int q_ = int(q / per);
      const int64_t rem = q % per;
      const int row = int(rem % Pc), col = int(rem / Pc);
      const int j = I.blocks[q_];
      const int r = block_row(row, K, Kj[j]);
      if(r >= 0)
        to_device(sdp.free_var_matrix.blocks[local_of[j]].LockedMatrix().CRef(r, col),
                  L, m, tot, q, e);
      else
        zero_device(L, m, tot, q, e);
    });
  }

  void upload_y(const SDP_Solver &solver)
  {
    // y: duplicated among the blocks
    std::vector<int32_t> m(size_t(L) * N);
    std::vector<int64_t> e(N);
    const auto &ym = solver.y.blocks.at(0).LockedMatrix();
    for(int i = 0; i < N; ++i)
      to_device(ym.CRef(i, 0), L, m.data(), N, i, e.data());
    gs->upload_y(m.data(), e.data());
  }

  void upload_x_X_Y(int ci, const SDP_Solver &solver)
  {
    const Chunk_Info &I = gs->chunk(ci);
    const int nb = I.nb, K = I.K, h = I.h, n = I.n, Pc = I.P;
    fill(ci, Gpu_Array::x, [&](int64_t q, int32_t *m, int64_t tot, int64_t *e) {
      const int q_ = int(q / Pc), row = int(q % Pc);
      const int j = I.blocks[q_];
      const int r = block_row(row, K, Kj[j]);
      if(r >= 0)
        to_device(solver.x.blocks[local_of[j]].LockedMatrix().CRef(r, 0), L, m, tot,
                  q, e);
      else
        zero_device(L, m, tot, q, e);
    });
    auto fill_matrix = [&](Gpu_Array which,
                           const Paired_Block_Diagonal_Matrix &M) {
      fill(ci, which, [&](int64_t q, int32_t *m, int64_t tot, int64_t *e) {
        const int64_t per = int64_t(n) * n;
        const int z = int(q / per);
        const int64_t rem = q % per;
        const int i0 = int(rem % n), j0 = int(rem / n);
        const int q_ = z % nb, par = z / nb;
        const int j = I.blocks[q_];
        const int hb = par == 0 ? he[j] : ho[j];
        const int ci_ = i0 / h, i = i0 % h, cj_ = j0 / h, jj = j0 % h;
        if(i < hb && jj < hb)
          to_device(M.get_block(local_of[j], par).LockedMatrix().CRef(
                      ci_ * hb + i, cj_ * hb + jj),
                    L, m, tot, q, e);
        else if(i0 == j0)
          one_device(L, m, tot, q, e);
        else
          zero_device(L, m, tot, q, e);
      });
    };
    fill_matrix(Gpu_Array::X, solver.X);
    fill_matrix(Gpu_Array::Y, solver.Y);
  }
};

Gpu_Compute::Gpu_Compute(const Solver_Parameters &parameters,
                         const Block_Info &block_info, const SDP &sdp,
                         const SDP_Solver &solver)
    : p_(std::make_unique<Impl>())
{
  ASSERT(El::mpi::Size(block_info.mpi_comm.value)==1,
         "--gpu requires complete blocks on each MPI rank");
  auto &P = *p_;
  P.L = limbs_for_precision(parameters.precision);
  P.N = int(sdp.dual_objective_b.Height());
  const int num_blocks = int(block_info.dimensions.size());
  std::vector<Block_Shape> shapes(num_blocks);
  P.local_of.assign(num_blocks, -1);
  for(size_t i = 0; i < block_info.block_indices.size(); ++i)
    P.local_of[block_info.block_indices[i]] = int(i);
  P.he.resize(num_blocks);
  P.ho.resize(num_blocks);
  P.Kj.resize(num_blocks);
  for(int j = 0; j < num_blocks; ++j)
    {
      const int loc = P.local_of[j];
      shapes[j].dim = int(block_info.dimensions[j]);
      shapes[j].K = int(block_info.num_points[j]);
      shapes[j].he = int(block_info.get_bilinear_bases_height(j,0));
      shapes[j].ho = int(block_info.get_bilinear_bases_height(j,1));
      P.he[j] = shapes[j].he;
      P.ho[j] = shapes[j].ho;
      P.Kj[j] = shapes[j].K;
      if(loc>=0)
        {
          ASSERT_EQUAL(sdp.bilinear_bases[2 * loc].Width(), shapes[j].K);
          ASSERT_EQUAL(sdp.bilinear_bases[2 * loc + 1].Width(), shapes[j].K);
          ASSERT_EQUAL(sdp.bilinear_bases[2 * loc].Height(), shapes[j].he);
          ASSERT_EQUAL(sdp.bilinear_bases[2 * loc + 1].Height(), shapes[j].ho);
        }
    }
  Gpu_Solver_Options o;
  o.world.rank=El::mpi::Rank();
  o.world.size=El::mpi::Size();
  o.world.device=-1;
  // the GPU by the rank's place on its node, as Block_Info's placement
  o.world.local_rank=block_info.node_comm.Rank();
  o.world.host_broadcast=[](void *p,size_t n,int root) {
    El::mpi::Broadcast(static_cast<unsigned char *>(p),int(n),root,El::mpi::COMM_WORLD);
  };
  o.q_tile=int(env_long("SDPB_GPU_Q_TILE",512));
  o.q_row_panel=int(env_long("SDPB_GPU_Q_ROW_PANEL",8192));
  o.q_wide_residual=env_long("SDPB_GPU_Q_WIDE_RESIDUAL",0)!=0;
  o.L = P.L;
  o.digits = int(env_long("SDPB_GPU_DIGITS", 4 * P.L));
  o.eig_tol = env_double("SDPB_GPU_EIG_TOL", 1e-10);
  o.correct_xy = env_long("SDPB_GPU_CORRECT_XY", 1) != 0;
  o.correct_schur = env_long("SDPB_GPU_CORRECT_SCHUR", 1) != 0;
  o.refine_dy = int(env_long("SDPB_GPU_REFINE_DY", 1));
  o.refine_system = int(env_long("SDPB_GPU_REFINE_SYSTEM", -1));
  o.refine_x = int(env_long("SDPB_GPU_REFINE_X", 0));
  o.schur_correction_digits = int(env_long("SDPB_GPU_SCHUR_DIGITS", -1));
  o.schur_correction_limit = int(env_long("SDPB_GPU_SCHUR_LIMIT", 0));
  o.crt_guard = int(env_long("SDPB_GPU_CRT_GUARD", 128));
  o.pq_guard = int(env_long("SDPB_GPU_PQ_GUARD", 0));
  o.structured_p = env_long("SDPB_GPU_STRUCTURED_P", 1) != 0;
  // 16: one chunk per dimension runs out of memory at nmax14 on one H200, with
  // 16 the peak is 131 GB. Not a planner yet.
  o.max_chunk_blocks = int(env_long("SDPB_GPU_CHUNK_BLOCKS", 16));
  o.max_chunk_rows = env_long("SDPB_GPU_CHUNK_ROWS", 0);
  o.sparse_b = env_long("SDPB_GPU_SPARSE_B", 0) != 0;
  o.lean_p = env_long("SDPB_GPU_LEAN_P", 0) != 0;
  o.product_timing = env_long("SDPB_GPU_PRODUCT_TIMING", 0) != 0;
  o.perturb_step = env_double("SDPB_GPU_PERTURB_STEP", 0);
  o.perturb_iteration = int(env_long("SDPB_GPU_PERTURB_ITERATION", 0));
  o.corrupt_leaf_limbs = int(env_long("SDPB_GPU_CORRUPT_LEAF_LIMBS", 0));
  o.corrupt_leaf_iteration = int(env_long("SDPB_GPU_CORRUPT_LEAF_ITERATION", 0));
  P.gs = std::make_unique<Gpu_Solver>(o, P.N, shapes);

  const auto t0 = std::chrono::steady_clock::now();
  const int L = P.L;
  P.upload_b(sdp);
  P.upload_y(solver);
  for(int ci = 0; ci < P.gs->num_chunks(); ++ci)
    {
      const Chunk_Info &I = P.gs->chunk(ci);
      for(int j:I.blocks) ASSERT(P.local_of[j]>=0,"GPU/host block ownership differs");
      const int nb = I.nb, K = I.K, h = I.h;
      // V (2 nb, h, K)
      P.fill(ci, Gpu_Array::V,
             [&](int64_t q, int32_t *m, int64_t tot, int64_t *e) {
               const int64_t per = int64_t(h) * K;
               const int z = int(q / per), rem = int(q % per);
               const int i = rem % h, k = rem / h;
               const int q_ = z % nb, par = z / nb;
               const auto &V = sdp.bilinear_bases[2 * P.local_of[I.blocks[q_]]
                                                  + par]
                                 .LockedMatrix();
               if(i < V.Height() && k < V.Width())
                 to_device(V.CRef(i, k), L, m, tot, q, e);
               else
                 zero_device(L, m, tot, q, e);
             });
      P.upload_c_B(ci, sdp);
      P.upload_x_X_Y(ci, solver);
    }
  P.gs->finish_load();
  if(El::mpi::Rank() == 0)
    {
      std::cout << "Gpu_Compute: " << P.gs->num_chunks() << " chunks, L = " << L
                << ", N = " << P.N << ", loaded in "
                << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - t0)
                     .count()
                << " s\n"
                << std::flush;
      for(int ci = 0; ci < P.gs->num_chunks(); ++ci)
        {
          const Chunk_Info &I = P.gs->chunk(ci);
          std::cout << "  chunk " << ci << ": dim " << I.dim << ", " << I.nb
                    << " blocks, K " << I.K << ", h " << I.h << ", P " << I.P
                    << ", n " << I.n << "\n";
        }
      std::cout << std::flush;
    }
}

Gpu_Compute::~Gpu_Compute() = default;

void Gpu_Compute::residues(const SDP &sdp, SDP_Solver &solver)
{
  auto &P = *p_;
  P.last = P.gs->residues();
  El::BigFloat psum(0);
  for(const auto &t : P.last.primal_terms)
    psum += from_host(t);
  if(El::mpi::Size()>1) psum=El::mpi::AllReduce(psum,El::mpi::SUM,El::mpi::COMM_WORLD);
  solver.primal_objective = sdp.objective_const + psum;
  solver.dual_objective = sdp.objective_const + from_host(P.last.dual_sum);
  solver.duality_gap
    = El::Abs(solver.primal_objective - solver.dual_objective)
      / El::Max(El::Abs(solver.primal_objective)
                  + El::Abs(solver.dual_objective),
                El::BigFloat(1));
  solver.primal_error_P = from_host(P.last.primal_error_P);
  solver.primal_error_p = from_host(P.last.primal_error_p);
  solver.dual_error = from_host(P.last.dual_error);
}

void Gpu_Compute::step(const Solver_Parameters &parameters, size_t total_psd_rows,
                       bool is_primal_and_dual_feasible, SDP_Solver &solver,
                       El::BigFloat &mu, El::BigFloat &beta_corrector,
                       El::BigFloat &primal_step_length,
                       El::BigFloat &dual_step_length, bool &terminate_now,
                       El::BigFloat &Q_cond_number,
                       El::BigFloat &max_block_cond_number,
                       std::string &max_block_cond_number_name)
{
  auto &P = *p_;
  const int L = P.L;
  El::BigFloat trace(0);
  for(const auto &t : P.last.minus_xy_trace)
    trace += from_host(t);
  if(El::mpi::Size()>1) trace=El::mpi::AllReduce(trace,El::mpi::SUM,El::mpi::COMM_WORLD);
  mu = -trace / total_psd_rows;
  if(mu > parameters.max_complementarity)
    {
      terminate_now = true;
      return;
    }
  P.gs->factor();
  solver.R_error = from_host(P.gs->R_error(to_host(mu, L)));

  const El::BigFloat beta_predictor
    = predictor_centering_parameter(parameters, is_primal_and_dual_feasible);
  P.gs->direction(to_host(beta_predictor * mu, L), false);

  // as corrector_centering_parameter() (step/centering_parameters.hxx)
  El::BigFloat frobenius(0);
  for(const auto &t : P.gs->frobenius_terms())
    frobenius += from_host(t);
  if(El::mpi::Size()>1) frobenius=El::mpi::AllReduce(frobenius,El::mpi::SUM,El::mpi::COMM_WORLD);
  const El::BigFloat r = frobenius / (mu * total_psd_rows);
  const El::BigFloat beta = r < 1 ? r * r : r;
  beta_corrector
    = is_primal_and_dual_feasible
        ? El::Min(El::Max(parameters.feasible_centering_parameter, beta),
                  El::BigFloat(1))
        : El::Max(parameters.infeasible_centering_parameter, beta);
  P.gs->direction(to_host(beta_corrector * mu, L), true);

  // update_cond_numbers.hxx, same order of blocks
  const Cond_Table &C = P.gs->cond();
  auto cond_of = [&](const Block_Cond &c) {
    if(c.max.m.empty() || c.min.m.empty())
      return El::BigFloat(0);
    const El::BigFloat ratio = from_host(c.max) / from_host(c.min);
    return ratio * ratio;
  };
  Q_cond_number = cond_of(C.Q);
  max_block_cond_number = 0;
  max_block_cond_number_name = "";
  int max_condition_order = std::numeric_limits<int>::max();
  for(size_t j = 0; j < C.schur.size(); ++j)
    {
      const El::BigFloat s = cond_of(C.schur[j]);
      if(max_block_cond_number < s)
        {
          max_block_cond_number = s;
          max_condition_order = int(5*j);
          max_block_cond_number_name
            = "schur_complement_cholesky.block_" + std::to_string(j);
        }
      for(int parity : {0, 1})
        {
          const El::BigFloat x = cond_of(C.X[2 * j + parity]);
          if(max_block_cond_number < x)
            {
              max_block_cond_number = x;
              max_condition_order = int(5*j+1+2*parity);
              max_block_cond_number_name = El::BuildString("X_cholesky.block_", j,
                                                           "_", parity);
            }
          const El::BigFloat y = cond_of(C.Y[2 * j + parity]);
          if(max_block_cond_number < y)
            {
              max_block_cond_number = y;
              max_condition_order = int(5*j+2+2*parity);
              max_block_cond_number_name = El::BuildString("Y_cholesky.block_", j,
                                                           "_", parity);
            }
        }
    }

  if(El::mpi::Size()>1)
    {
      const auto local=max_block_cond_number;
      max_block_cond_number=El::mpi::AllReduce(local,El::mpi::MAX,El::mpi::COMM_WORLD);
      // Preserve the single-rank block traversal on exact ties, regardless
      // of which rank owns the earliest block in that traversal.
      const int order=local==max_block_cond_number ? max_condition_order
        : std::numeric_limits<int>::max();
      const int first=El::mpi::AllReduce(order,El::mpi::MIN,El::mpi::COMM_WORLD);
      int root=order==first ? El::mpi::Rank() : El::mpi::Size();
      root=El::mpi::AllReduce(root,El::mpi::MIN,El::mpi::COMM_WORLD);
      int length=int(max_block_cond_number_name.size());
      El::mpi::Broadcast(&length,1,root,El::mpi::COMM_WORLD);
      max_block_cond_number_name.resize(length);
      El::mpi::Broadcast(reinterpret_cast<unsigned char *>(max_block_cond_number_name.data()),
                         length,root,El::mpi::COMM_WORLD);
    }

  // a block that the float64 eigenvalue cannot resolve: SDPB's multiprecision
  // eigenvalue (min_eigenvalue(), step/step_length.hxx) on the host
  auto refine = [&](int n, const std::vector<int32_t> &m,
                    const std::vector<int64_t> &e) -> Host_Float {
    El::Matrix<El::BigFloat> A(n, n);
    const int64_t per = int64_t(n) * n;
    for(int j = 0; j < n; ++j)
      for(int i = 0; i < n; ++i)
        {
          const int64_t q = int64_t(j) * n + i;
          A.Set(i, j,
                e[q] == MPF_ZERO_EXP
                  ? El::BigFloat(0)
                  : from_limbs(reinterpret_cast<const uint32_t *>(m.data()) + q,
                               L, per, e[q]));
        }
    El::Matrix<El::BigFloat> w;
    El::HermitianEigCtrl<El::BigFloat> ctrl;
    ctrl.tridiagEigCtrl.dcCtrl.cutoff = n / 2 + 1;
    ctrl.tridiagEigCtrl.dcCtrl.secularCtrl.maxIterations = 16384;
    El::HermitianEig(El::UpperOrLowerNS::LOWER, A, w, ctrl);
    El::BigFloat lambda = w.CRef(0, 0);
    for(int i = 1; i < w.Height(); ++i)
      lambda = El::Min(lambda, w.CRef(i, 0));
    return to_host(lambda, L);
  };

  const double gamma = double(parameters.step_length_reduction);
  const Min_Eig lamP = P.gs->min_eig(false, gamma, refine);
  const Min_Eig lamD = P.gs->min_eig(true, gamma, refine);
  auto step_of = [&](const Min_Eig &m) {
    if(!m.any)
      return El::BigFloat(1);
    const El::BigFloat lam = from_host(m.value);
    if(lam > -parameters.step_length_reduction)
      return El::BigFloat(1);
    return -parameters.step_length_reduction / lam;
  };
  primal_step_length = step_of(lamP);
  dual_step_length = step_of(lamD);
  if(El::mpi::Size()>1)
    {
      primal_step_length=El::mpi::AllReduce(primal_step_length,El::mpi::MIN,El::mpi::COMM_WORLD);
      dual_step_length=El::mpi::AllReduce(dual_step_length,El::mpi::MIN,El::mpi::COMM_WORLD);
    }
  if(lamP.refined + lamD.refined > 0)
    std::cerr << "Gpu_Compute: " << lamP.refined + lamD.refined
              << " blocks' step eigenvalues refined in multiprecision\n";
  if(is_primal_and_dual_feasible)
    {
      primal_step_length = El::Min(primal_step_length, dual_step_length);
      dual_step_length = primal_step_length;
    }
  P.gs->update(to_host(primal_step_length, L), to_host(dual_step_length, L));
  P.gs->end_iteration();
}

void Gpu_Compute::get_y(SDP_Solver &solver)
{
  auto &P = *p_;
  std::vector<int32_t> m(size_t(P.L) * P.N);
  std::vector<int64_t> e(P.N);
  P.gs->download_y(m.data(), e.data());
  for(auto &yb : solver.y.blocks)
    {
      auto &Y = yb.Matrix();
      for(int i = 0; i < P.N; ++i)
        Y.Set(i, 0,
              e[i] == MPF_ZERO_EXP
                ? El::BigFloat(0)
                : from_limbs(reinterpret_cast<const uint32_t *>(m.data()) + i,
                             P.L, P.N, e[i]));
    }
}

void Gpu_Compute::set_iterate(const SDP_Solver &solver)
{
  auto &P = *p_;
  P.upload_y(solver);
  for(int ci = 0; ci < P.gs->num_chunks(); ++ci)
    P.upload_x_X_Y(ci, solver);
}

void Gpu_Compute::end_run() { p_->gs->end_run(); }

void Gpu_Compute::update_blocks(const SDP &sdp, const std::vector<size_t> &changed,
                                const bool objectives_changed)
{
  auto &P = *p_;
  if(objectives_changed)
    P.upload_b(sdp);
  for(int ci = 0; ci < P.gs->num_chunks(); ++ci)
    {
      const auto &blocks = P.gs->chunk(ci).blocks;
      const bool touched
        = std::any_of(blocks.begin(), blocks.end(), [&](const int j) {
            return std::find(changed.begin(), changed.end(), size_t(j))
                   != changed.end();
          });
      if(touched)
        P.upload_c_B(ci, sdp);
    }
}

void Gpu_Compute::get_iterate(SDP_Solver &solver)
{
  auto &P = *p_;
  const int L = P.L;
  auto get = [&](int ci, Gpu_Array a, std::vector<int32_t> &m,
                 std::vector<int64_t> &e) {
    const int64_t tot = P.gs->entries(ci, a);
    m.resize(size_t(L) * tot);
    e.resize(tot);
    P.gs->download(ci, a, m.data(), e.data());
    return tot;
  };
  auto bf = [&](const std::vector<int32_t> &m, const std::vector<int64_t> &e,
                int64_t tot, int64_t q) {
    return e[q] == MPF_ZERO_EXP
             ? El::BigFloat(0)
             : from_limbs(reinterpret_cast<const uint32_t *>(m.data()) + q, L,
                          tot, e[q]);
  };
  {
    std::vector<int32_t> m(size_t(L) * P.N);
    std::vector<int64_t> e(P.N);
    P.gs->download_y(m.data(), e.data());
    for(auto &yb : solver.y.blocks)
      {
        auto &Y = yb.Matrix();
        for(int i = 0; i < P.N; ++i)
          Y.Set(i, 0, bf(m, e, P.N, i));
      }
  }
  std::vector<int32_t> m;
  std::vector<int64_t> e;
  for(int ci = 0; ci < P.gs->num_chunks(); ++ci)
    {
      const Chunk_Info &I = P.gs->chunk(ci);
      const int nb = I.nb, K = I.K, h = I.h, n = I.n, Pc = I.P;
      int64_t tot = get(ci, Gpu_Array::x, m, e);
      for(int q_ = 0; q_ < nb; ++q_)
        {
          const int j = I.blocks[q_];
          auto &xb = solver.x.blocks[P.local_of[j]].Matrix();
          const int Kb = P.Kj[j];
          const int nseg = Pc / K;
          for(int s = 0; s < nseg; ++s)
            for(int k = 0; k < Kb; ++k)
              xb.Set(s * Kb + k, 0,
                     bf(m, e, tot, int64_t(q_) * Pc + s * K + k));
        }
      auto get_matrix = [&](Gpu_Array which,
                            Paired_Block_Diagonal_Matrix &M) {
        const int64_t t2 = get(ci, which, m, e);
        const int64_t per = int64_t(n) * n;
        for(int z = 0; z < 2 * nb; ++z)
          {
            const int q_ = z % nb, par = z / nb;
            const int j = I.blocks[q_];
            const int hb = par == 0 ? P.he[j] : P.ho[j];
            auto &Mb = M.get_block(P.local_of[j], par).Matrix();
            parallel_for(per, [&](int64_t q) {
              const int i0 = int(q % n), j0 = int(q / n);
              const int i = i0 % h, jj = j0 % h;
              if(i < hb && jj < hb)
                Mb.Set((i0 / h) * hb + i, (j0 / h) * hb + jj,
                       bf(m, e, t2, z * per + q));
            });
          }
      };
      get_matrix(Gpu_Array::X, solver.X);
      get_matrix(Gpu_Array::Y, solver.Y);
    }
}

void Gpu_Compute::write_profile(const std::filesystem::path &path) const
{
  const auto dest=El::mpi::Rank()==0 ? path : std::filesystem::path(
    path.string()+".rank"+std::to_string(El::mpi::Rank()));
  std::ofstream os(dest);
  os << p_->gs->profile_json();
}
