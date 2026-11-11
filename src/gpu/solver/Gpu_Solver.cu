// GPU backend for SDPB: see Gpu_Solver.hxx.  Based on SDPviBe by Antoine Tilloy
// (https://github.com/AntoineTilloy/SDPviBe), the prototype of this backend:
// the passes below follow sdpvibe/solver.py (Solver.run, cmws, pairings, schur,
// seg_diag, min_eig), with explicit inverses in place of triangular solves
// (linalg/Block_Factor).  Portions adapted from SDPviBe, Copyright (c) 2026
// Antoine Tilloy, MIT License (docs/GPU_backend_credits.md).

#include "Gpu_Solver.hxx"
#include "Gpu_Layout.hxx"
#include "../linalg/Distributed_Q.hxx"
#include "Mpf_Ops.hxx"
#include "../rns/rns_kernels.cuh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

#include <sys/resource.h>

#define GS_CHECK(x)                                                           \
  do                                                                          \
    {                                                                         \
      cudaError_t e_ = (x);                                                   \
      if(e_ != cudaSuccess)                                                   \
        {                                                                     \
          std::ostringstream ss_;                                             \
          ss_ << __FILE__ << ":" << __LINE__ << ": "                          \
              << cudaGetErrorString(e_);                                      \
          throw std::runtime_error(ss_.str());                                \
        }                                                                     \
  } while(0)

extern "C" void dsyevd_(const char *jobz, const char *uplo, const int *n,
                        double *a, const int *lda, double *w, double *work,
                        const int *lwork, int *iwork, const int *liwork,
                        int *info);

namespace sdpb_gpu
{
  namespace
  {
    typedef long double ld;

    // the products' work buffers, shared by all engines: one product runs at a
    // time, so the largest product sets them, not the sum over shapes and batch
    // sizes (Rns_Workspace)
    Rns_Workspace *g_workspace = nullptr;

    cudaError_t solver_malloc(void **p, size_t bytes)
    {
      cudaError_t e = malloc_async_retry(p, bytes, 0);
      if(e == cudaErrorMemoryAllocation && g_workspace)
        {
          (void)cudaGetLastError();
          cudaStreamSynchronize(0);
          g_workspace->release();
          e = malloc_async_retry(p, bytes, 0);
        }
      return e;
    }

    // An owning device array
    struct Mat
    {
      Mpf_Batch b;
      Mat() = default;
      Mat(int L, int batch, int64_t rows, int64_t cols)
      {
        b.L = L;
        b.batch = batch;
        b.rows = rows;
        b.cols = cols;
        GS_CHECK(solver_malloc(reinterpret_cast<void **>(&b.m),
                               size_t(L) * b.entries() * 4));
        GS_CHECK(solver_malloc(reinterpret_cast<void **>(&b.e),
                               size_t(b.entries()) * 8));
      }
      Mat(const Mat &) = delete;
      Mat &operator=(const Mat &) = delete;
      Mat(Mat &&o) noexcept : b(o.b) { o.b = Mpf_Batch{}; }
      Mat &operator=(Mat &&o) noexcept
      {
        if(this != &o)
          {
            free();
            b = o.b;
            o.b = Mpf_Batch{};
          }
        return *this;
      }
      ~Mat() { free(); }
      void free()
      {
        if(b.m)
          cudaFreeAsync(b.m, 0);
        if(b.e)
          cudaFreeAsync(b.e, 0);
        b.m = nullptr;
        b.e = nullptr;
      }
      explicit operator bool() const { return b.m != nullptr; }
    };

    struct Chunk
    {
      Chunk_Info info;
      int nseg = 0;
      std::vector<int> hj, Kj;
      int *d_hj = nullptr, *d_Kj = nullptr;
      Mat V, VT, VTrep, Vx, c, B, x, X, Y;
      Mat Bnz;                 // with sparse_b: B's nonzero entries; B only during a use
      uint32_t *Bpos = nullptr; // their flat indices in B
      Mat Xi, XiT, Yi, YiT, AXI, AY, dres, Pres, mXY;
      Mat Si, SiT, Pf;
      Mat R, w, dx, dX, dY, S, rhs, w2;
    };

    Mat make_copy(const Mat &A)
    {
      Mat C(A.b.L, A.b.batch, A.b.rows, A.b.cols);
      GS_CHECK(cudaMemcpyAsync(C.b.m, A.b.m, size_t(A.b.L) * A.b.entries() * 4,
                               cudaMemcpyDeviceToDevice, 0));
      GS_CHECK(cudaMemcpyAsync(C.b.e, A.b.e, size_t(A.b.entries()) * 8,
                               cudaMemcpyDeviceToDevice, 0));
      return C;
    }

    Mat make_zero(int L, int batch, int64_t rows, int64_t cols)
    {
      Mat C(L, batch, rows, cols);
      mpf_zero(C.b, 0);
      return C;
    }

    bool mag_less(const uint32_t *am, int64_t ae, const uint32_t *bm, int64_t be,
                  int L)
    {
      if(ae != be)
        return ae < be;
      for(int t = L - 1; t >= 0; --t)
        if(am[t] != bm[t])
          return am[t] < bm[t];
      return false;
    }
  }

  struct Gpu_Solver::Impl
  {
    Gpu_Solver_Options opt;
    int N, L, D;
    std::vector<Block_Shape> shapes;
    std::vector<Chunk> chunks;
    int nblocks_total;
    Gram_Options gopt;
    std::unique_ptr<Rns_Workspace> workspace;
    std::map<int, std::unique_ptr<Block_Factor>> engines;
    std::unique_ptr<Mpf_Gram> gram;
    std::unique_ptr<Gpu_Collective> world;
    std::unique_ptr<Distributed_Q> distributed_q;
    int64_t global_rows = 0;
    std::vector<std::string> distributed_profiles;
    Mat b, y, dy, pvec, Xq, Qfull;
    Cond_Table conds;
    int iteration = 1;
    int refine_system_now = 0;

    // profile: seconds by stage, this iteration and in total
    using Clock = std::chrono::steady_clock;
    Clock::time_point t_stage;
    std::string cur_stage;
    std::map<std::string, double> total_s, iter_s, last_iter_s, stage_used_gb,
      stage_host_gb;
    std::map<std::string, uint64_t> total_payload, iter_payload, last_iter_payload;
    std::map<std::string, double> total_comm_s, iter_comm_s, last_iter_comm_s;
    double stage_comm_begin = 0;
    uint64_t stage_payload_begin = 0;
    std::vector<double> iter_wall;
    Clock::time_point t_iter;
    double peak_bytes = 0;

    // High-water marks of the stage just ended: the stream-ordered pool's
    // used-memory high (not its reserve, which nvidia-smi would show), plus what
    // cudaMalloc holds outside the pool; and the host's peak resident size so far.
    void note_memory()
    {
      int dev = 0;
      cudaMemPool_t pool;
      GS_CHECK(cudaGetDevice(&dev));
      GS_CHECK(cudaDeviceGetDefaultMemPool(&pool, dev));
      size_t free_b = 0, total_b = 0;
      GS_CHECK(cudaMemGetInfo(&free_b, &total_b));
      unsigned long long used_high = 0, reserved = 0;
      GS_CHECK(cudaMemPoolGetAttribute(pool, cudaMemPoolAttrUsedMemHigh, &used_high));
      GS_CHECK(cudaMemPoolGetAttribute(pool, cudaMemPoolAttrReservedMemCurrent,
                                       &reserved));
      const double outside_pool
        = std::max(0.0, double(total_b - free_b) - double(reserved));
      const double high = double(used_high) + outside_pool;
      peak_bytes = std::max(peak_bytes, high);
      stage_used_gb[cur_stage] = std::max(stage_used_gb[cur_stage], high / 1e9);
      unsigned long long zero = 0;
      GS_CHECK(cudaMemPoolSetAttribute(pool, cudaMemPoolAttrUsedMemHigh, &zero));
      struct rusage ru;
      getrusage(RUSAGE_SELF, &ru);
      stage_host_gb[cur_stage]
        = std::max(stage_host_gb[cur_stage], double(ru.ru_maxrss) * 1024.0 / 1e9);
    }

    void stage(const char *name)
    {
      GS_CHECK(cudaStreamSynchronize(0));
      const auto now = Clock::now();
      if(!cur_stage.empty())
        {
          const double dt = std::chrono::duration<double>(now - t_stage).count();
          total_s[cur_stage] += dt;
          iter_s[cur_stage] += dt;
          if(world)
            {
              const uint64_t bytes=world->received_bytes()-stage_payload_begin;
              total_payload[cur_stage]+=bytes;
              iter_payload[cur_stage]+=bytes;
              const double comm=world->seconds()-stage_comm_begin;
              total_comm_s[cur_stage]+=comm;
              iter_comm_s[cur_stage]+=comm;
            }
          note_memory();
        }
      cur_stage = name ? name : "";
      stage_payload_begin=world ? world->received_bytes() : 0;
      stage_comm_begin=world ? world->seconds() : 0;
      t_stage = Clock::now();
    }

    Block_Factor &eng(int batch)
    {
      auto it = engines.find(batch);
      if(it == engines.end())
        it = engines
               .emplace(batch, std::make_unique<Block_Factor>(32, batch, L, D,
                                                              gopt,
                                                              cudaStream_t(0)))
               .first;
      return *it->second;
    }

    // C = A^T B
    Mat prod(const Mat &A, const Mat &B)
    {
      Mat C(L, A.b.batch, A.b.cols, B.b.cols);
      eng(A.b.batch).mm(A.b, true, B.b, false, C.b);
      return C;
    }

    Mat add(const Mat &A, const Mat &B, int sgn)
    {
      Mat C(L, A.b.batch, A.b.rows, A.b.cols);
      mpf_batch_add(A.b, B.b, sgn, C.b, 0);
      return C;
    }

    Mat transpose(const Mat &A)
    {
      Mat T(L, A.b.batch, A.b.cols, A.b.rows);
      mpf_transpose(A.b, T.b, 0);
      return T;
    }

    // (A + A^T) / 2
    Mat symmetrize(const Mat &A)
    {
      Mat T = transpose(A), S(L, A.b.batch, A.b.rows, A.b.cols);
      mpf_batch_add(A.b, T.b, 1, S.b, 0);
      mpf_shift_exp(S.b, -1, 0);
      return S;
    }

    Mat scalar(const Host_Float &v)
    {
      Mat S(L, 1, 1, 1);
      std::vector<int32_t> m(L, 0);
      int64_t e = -(int64_t(1) << 40);
      const size_t n = v.m.size();
      if(n)
        {
          // normalize to exactly NB bits (the host values already are, except
          // sums with more limbs)
          const bool neg = int32_t(v.m[n - 1]) < 0;
          std::vector<uint32_t> a(v.m);
          if(neg)
            {
              unsigned long long c = 1;
              for(auto &w : a)
                {
                  const unsigned long long t = (unsigned long long)(~w) + c;
                  w = uint32_t(t);
                  c = t >> 32;
                }
            }
          int64_t bl = 0;
          for(int k = int(n) - 1; k >= 0; --k)
            if(a[k])
              {
                bl = 32 * k + 32 - __builtin_clz(a[k]);
                break;
              }
          if(bl)
            {
              const int64_t NB = 32 * L - 2;
              int64_t sh = bl - NB;
              auto extract = [&](int64_t s, std::vector<uint32_t> &r) {
                r.assign(L + 1, 0);
                for(int t = 0; t <= L; ++t)
                  {
                    const int64_t bit = 32 * t + s;
                    const int64_t q = bit >> 5;
                    const int rb = int(bit & 31);
                    const uint32_t lo = (q >= 0 && q < int64_t(n)) ? a[q] : 0u;
                    const uint32_t hi
                      = (q + 1 >= 0 && q + 1 < int64_t(n)) ? a[q + 1] : 0u;
                    r[t] = uint32_t((((unsigned long long)hi << 32) | lo) >> rb);
                  }
                if(s >= 1)
                  {
                    const int64_t gb = s - 1;
                    unsigned long long c
                      = (gb >> 5) < int64_t(n)
                          ? (a[gb >> 5] >> (gb & 31)) & 1u
                          : 0u;
                    for(int t = 0; t <= L && c; ++t)
                      {
                        const unsigned long long w
                          = (unsigned long long)r[t] + c;
                        r[t] = uint32_t(w);
                        c = w >> 32;
                      }
                  }
              };
              std::vector<uint32_t> r;
              extract(sh, r);
              if(r[L] != 0 || (r[L - 1] >> (NB - 32 * (L - 1))) != 0)
                {
                  ++sh;
                  extract(sh, r);
                }
              if(neg)
                {
                  unsigned long long c = 1;
                  for(int t = 0; t < L; ++t)
                    {
                      const unsigned long long w
                        = (unsigned long long)(~r[t]) + c;
                      r[t] = uint32_t(w);
                      c = w >> 32;
                    }
                }
              for(int t = 0; t < L; ++t)
                m[t] = int32_t(r[t]);
              e = v.e + sh;
            }
        }
      GS_CHECK(cudaStreamSynchronize(0));
      GS_CHECK(cudaMemcpy(S.b.m, m.data(), 4 * size_t(L), cudaMemcpyHostToDevice));
      GS_CHECK(cudaMemcpy(S.b.e, &e, 8, cudaMemcpyHostToDevice));
      return S;
    }

    Host_Float max_of(const Host_Float &a, const Host_Float &b) const
    {
      // both are nonnegative magnitudes rounded to floats
      if(a.m.empty() || a.e == -(int64_t(1) << 40))
        return b;
      if(b.m.empty() || b.e == -(int64_t(1) << 40))
        return a;
      return mag_less(a.m.data(), a.e, b.m.data(), b.e, L) ? b : a;
    }

    Host_Float global_max(const Host_Float &value) const
    {
      if(!world) return value;
      Host_Float result;
      for(int rank=0;rank<world->size();++rank)
        {
          Host_Float v=value;
          v.m.resize(L,0);
          if(value.m.empty()) v.e=-(int64_t(1)<<40);
          opt.world.host_broadcast(v.m.data(),size_t(L)*4,rank);
          opt.world.host_broadcast(&v.e,8,rank);
          result=max_of(result,v);
        }
      return result;
    }
    Mat root_copy(const Mat &v)
    {
      return !world || world->rank()==0 ? make_copy(v)
        : make_zero(L,v.b.batch,v.b.rows,v.b.cols);
    }
    void sum_ranks(Mat &v)
    {
      if(!world) return;
      Mat sum(L,v.b.batch,v.b.rows,v.b.cols);
      distributed_sum(v.b,sum.b,*world);
      v=std::move(sum);
    }
    bool want_Q_digest() const
    { return iteration==1 && std::getenv("SDPB_GPU_Q_DIGEST_PREFIX"); }
    void record_Q_digest(const std::string &digest)
    {
      const char *prefix=std::getenv("SDPB_GPU_Q_DIGEST_PREFIX");
      std::ofstream out(std::string(prefix)+".rank"+std::to_string(opt.world.rank)+".txt");
      out << digest;
      if(!out) throw std::runtime_error("could not write Q tile digests");
    }

    void record_distributed_Q_peak()
    {
      const double high=double(distributed_q->last_operation_peak_device_bytes());
      peak_bytes=std::max(peak_bytes,high);
      stage_used_gb[cur_stage]=std::max(stage_used_gb[cur_stage],high/1e9);
    }

    void distributed_Q_factor()
    {
      stage("Q = P^T P");
      Gram_Options qopt=gopt;
      qopt.timing=opt.product_timing;
      distributed_q=std::make_unique<Distributed_Q>(N,L,D,global_rows,*world,qopt,
        Distributed_Q_Options{opt.q_tile,opt.q_row_panel,opt.refine_dy,opt.q_wide_residual});
      std::vector<Mpf_Matrix> panels;
      for(auto &c:chunks) panels.push_back(P_of(c).view());
      distributed_q->assemble(panels);
      record_distributed_Q_peak();
      if(opt.lean_p) for(auto &c:chunks) c.Pf.free();
      if(want_Q_digest()) record_Q_digest(distributed_q->digests());
      stage("Cholesky Q + inverse");
      distributed_q->factor();
      record_distributed_Q_peak();
      Mat diag(L,1,N,1);
      distributed_q->factor_diagonal(diag.b);
      std::vector<int32_t> m(size_t(L)*N);
      std::vector<int64_t> e(N);
      mpf_download(diag.b,m.data(),e.data(),0);
      Block_Cond c;
      for(int i=0;i<N;++i)
        {
          Host_Float v; v.e=e[i]; v.m.resize(L);
          for(int t=0;t<L;++t) v.m[t]=uint32_t(m[size_t(t)*N+i]);
          if(e[i]==-(int64_t(1)<<40) || int32_t(v.m.back())<=0)
            throw std::runtime_error("distributed Q has a nonpositive factor diagonal");
          c.max=max_of(c.max,v);
          if(c.min.m.empty() || mag_less(v.m.data(),v.e,c.min.m.data(),c.min.e,L)) c.min=v;
        }
      conds.Q=c;
      peak_bytes=std::max(peak_bytes,double(distributed_q->peak_device_bytes()));
    }

    // (L_S, X) with X = L_S^-1 of the padded A: returns the inverse cut to n x n
    // and its transpose, and records the diagonal of the factor.  valid(z, i)
    // says which rows of block z count in the condition number.
    template <class Valid>
    void factor_block(const Mat &A, Mat &Xi, Mat *XiT, Factor_Options fo,
                      std::vector<Block_Cond> &out,
                      const std::vector<size_t> &out_index, const Valid &valid)
    {
      const int n = int(A.b.rows), batch = A.b.batch;
      const int np = std::max(32, (n + 31) / 32 * 32);
      Block_Factor &E = eng(batch);
      Mat Ap;
      const Mpf_Batch *in = &A.b;
      if(np != n)
        {
          Ap = Mat(L, batch, np, np);
          mpf_pad_identity(A.b, Ap.b, 0);
          in = &Ap.b;
        }
      Mat Lf(L, batch, np, np), Xp(L, batch, np, np);
      E.factor(*in, Lf.b, Xp.b, fo);
      Mat D(L, batch, np, 1);
      mpf_diagonal(Lf.b, D.b, 0);
      std::vector<int32_t> dm(size_t(L) * D.b.entries());
      std::vector<int64_t> de(D.b.entries());
      mpf_download(D.b, dm.data(), de.data(), 0);
      const int64_t tot = D.b.entries();
      for(int z = 0; z < batch; ++z)
        {
          int64_t imax = -1, imin = -1;
          std::vector<uint32_t> cm(L), bmax(L), bmin(L);
          for(int i = 0; i < np; ++i)
            {
              if(!valid(z, i))
                continue;
              const int64_t q = int64_t(z) * np + i;
              for(int t = 0; t < L; ++t)
                cm[t] = uint32_t(dm[int64_t(t) * tot + q]);
              if(de[q] == -(int64_t(1) << 40) || int32_t(cm[L - 1]) <= 0)
                throw std::runtime_error(
                  "Gpu_Solver: Cholesky factor has a nonpositive diagonal "
                  "entry: matrix is not positive definite");
              if(imax < 0 || mag_less(bmax.data(), de[imax], cm.data(), de[q], L))
                {
                  imax = q;
                  bmax = cm;
                }
              if(imin < 0 || mag_less(cm.data(), de[q], bmin.data(), de[imin], L))
                {
                  imin = q;
                  bmin = cm;
                }
            }
          Block_Cond c;
          if(imax >= 0)
            {
              c.max.m = bmax;
              c.max.e = de[imax];
              c.min.m = bmin;
              c.min.e = de[imin];
            }
          out[out_index[z]] = c;
        }
      if(np != n)
        {
          Xi = Mat(L, batch, n, n);
          mpf_sub_block(Xp.b, 0, 0, Xi.b, false, 0);
        }
      else
        Xi = std::move(Xp);
      if(XiT)
        *XiT = transpose(Xi);
    }

    // sum_p a_p A_p for a (nb, P, 1): R (2 nb, n, n)
    Mat cmws(Chunk &c, const Mat &a)
    {
      const auto &I = c.info;
      Mat VD(L, 2 * I.nb * c.nseg, I.K, I.h);
      mpf_scale_segments(c.VT.b, a.b, VD.b, c.nseg, 0);
      Mat C = prod(c.VTrep, VD);
      VD.free();
      Mat R(L, 2 * I.nb, I.n, I.n);
      mpf_cmws_assemble(C.b, R.b, I.h, I.dim, 0);
      return R;
    }

    // the one place that says where a chunk's B and P are (streaming later
    // changes only these)
    const Mpf_Batch &B_of(Chunk &c)
    {
      if(opt.sparse_b && !c.B)
        {
          c.B = Mat(L, c.info.nb, c.info.P, N);
          mpf_unpack_sparse(c.Bnz.b, c.Bpos, c.B.b, 0);
        }
      return c.B.b;
    }
    // after each use of B_of: only a sparse chunk's expansion goes
    void B_done(Chunk &c)
    {
      if(opt.sparse_b)
        c.B.free();
    }
    void upload_sparse_B(Chunk &c, const int32_t *m, const int64_t *e)
    {
      std::vector<uint32_t> pos;
      std::vector<int32_t> vm;
      std::vector<int64_t> ve;
      mpf_pack_sparse(m, e, c.B.b.entries(), L, pos, vm, ve);
      // the chunk's earlier data, if this replaces them (resident mode)
      c.Bnz.free();
      if(c.Bpos)
        {
          GS_CHECK(cudaFreeAsync(c.Bpos, 0));
          c.Bpos = nullptr;
        }
      c.Bnz = Mat(L, 1, int64_t(pos.size()), 1);
      mpf_upload(vm.data(), ve.data(), c.Bnz.b, 0);
      GS_CHECK(solver_malloc(reinterpret_cast<void **>(&c.Bpos),
                             pos.size() * sizeof(uint32_t)));
      GS_CHECK(cudaMemcpy(c.Bpos, pos.data(), pos.size() * sizeof(uint32_t),
                          cudaMemcpyHostToDevice));
    }
    static const Mpf_Batch &P_of(const Chunk &c) { return c.Pf.b; }

    Mat gemv_B_y(Chunk &c, const Mat &v)
    {
      // B v: (nb, P, 1)
      const auto &I = c.info;
      Gemv_Spec g;
      g.batch = I.nb;
      g.no = I.P;
      g.ab = I.P * int64_t(N);
      g.ao = 1;
      g.ati = I.P;
      g.nti = N;
      g.vti = 1;
      g.ob = I.P;
      g.oo = 1;
      Mat O(L, I.nb, I.P, 1);
      mpf_gemv(B_of(c), v.b, O.b, g, 0);
      B_done(c);
      return O;
    }

    Mat gemv_Bt(Chunk &c, const Mat &v)
    {
      // B^T v: (1, N, 1), v (nb, P, 1)
      const auto &I = c.info;
      Gemv_Spec g;
      g.batch = 1;
      g.no = N;
      g.ao = I.P;
      g.atb = I.P * int64_t(N);
      g.ati = 1;
      g.ntb = I.nb;
      g.nti = I.P;
      g.vtb = I.P;
      g.vti = 1;
      g.ob = 0;
      g.oo = 1;
      Mat O(L, 1, N, 1);
      mpf_gemv(B_of(c), v.b, O.b, g, 0);
      B_done(c);
      return O;
    }

    Mat gemv_Pt(Chunk &c, const Mat &v)
    {
      // P^T v, P flat (nb P, N): (1, N, 1)
      if(opt.lean_p)
        return gemv_Bt(c, gemv_square(c.Si, v, true)); // B^T L^-T v
      const int64_t R = int64_t(c.info.nb) * c.info.P;
      Gemv_Spec g;
      g.no = N;
      g.ao = R;
      g.ati = 1;
      g.nti = R;
      g.vti = 1;
      g.oo = 1;
      Mat O(L, 1, N, 1);
      mpf_gemv(P_of(c), v.b, O.b, g, 0);
      return O;
    }

    Mat gemv_P(Chunk &c, const Mat &v)
    {
      // P v: (nb, P, 1), the flat rows are the batch layout of a vector
      if(opt.lean_p)
        return gemv_square(c.Si, gemv_B_y(c, v), false); // L^-1 B v
      const int64_t R = int64_t(c.info.nb) * c.info.P;
      Gemv_Spec g;
      g.no = R;
      g.ao = 1;
      g.ati = R;
      g.nti = N;
      g.vti = 1;
      g.oo = 1;
      Mat O(L, c.info.nb, c.info.P, 1);
      mpf_gemv(P_of(c), v.b, O.b, g, 0);
      return O;
    }

    // M v per block, M (batch, P, P), v (batch, P, 1); transposed: M^T v
    Mat gemv_square(const Mat &M, const Mat &v, bool transposed)
    {
      const int64_t P = M.b.rows;
      Gemv_Spec g;
      g.batch = M.b.batch;
      g.no = P;
      g.ab = P * P;
      g.ao = transposed ? P : 1;
      g.ati = transposed ? 1 : P;
      g.nti = P;
      g.vb = P;
      g.vti = 1;
      g.ob = P;
      g.oo = 1;
      Mat O(L, M.b.batch, P, 1);
      mpf_gemv(M.b, v.b, O.b, g, 0);
      return O;
    }

    void build_chunks()
    {
      const auto layout=gpu_chunks(shapes,opt.max_chunk_blocks,opt.max_chunk_rows);
      const auto owners=gpu_chunk_owners(layout,opt.world.size);
      global_rows=0;
      for(size_t ci=0;ci<layout.size();++ci)
        {
          global_rows+=int64_t(layout[ci].nb)*layout[ci].P;
          if(owners[ci]!=opt.world.rank) continue;
          Chunk c; c.info=layout[ci];
          auto &I=c.info;
          c.nseg=I.dim*(I.dim+1)/2;
          c.hj.resize(2 * I.nb);
          c.Kj.resize(I.nb);
          for(int q = 0; q < I.nb; ++q)
            {
              const auto &s = shapes[I.blocks[q]];
              c.hj[q] = s.he;
              c.hj[I.nb + q] = s.ho;
              c.Kj[q] = s.K;
            }
          GS_CHECK(cudaMalloc(&c.d_hj, c.hj.size() * sizeof(int)));
          GS_CHECK(cudaMalloc(&c.d_Kj, c.Kj.size() * sizeof(int)));
          GS_CHECK(cudaMemcpy(c.d_hj, c.hj.data(), c.hj.size() * sizeof(int),
                              cudaMemcpyHostToDevice));
          GS_CHECK(cudaMemcpy(c.d_Kj, c.Kj.data(), c.Kj.size() * sizeof(int),
                              cudaMemcpyHostToDevice));
          c.V = Mat(L, 2 * I.nb, I.h, I.K);
          c.c = Mat(L, I.nb, I.P, 1);
          if(opt.sparse_b)
            {
              c.B.b.L = L; // the shape only; upload keeps the nonzeros
              c.B.b.batch = I.nb;
              c.B.b.rows = I.P;
              c.B.b.cols = N;
            }
          else
            c.B = Mat(L, I.nb, I.P, N);
          c.x = Mat(L, I.nb, I.P, 1);
          c.X = Mat(L, 2 * I.nb, I.n, I.n);
          c.Y = Mat(L, 2 * I.nb, I.n, I.n);
          chunks.push_back(std::move(c));
        }
      nblocks_total = int(shapes.size());
      conds.schur.resize(nblocks_total);
      conds.X.resize(2 * nblocks_total);
      conds.Y.resize(2 * nblocks_total);
    }

    Mat &array(int i, Gpu_Array a)
    {
      auto &c = chunks.at(i);
      switch(a)
        {
        case Gpu_Array::V: return c.V;
        case Gpu_Array::c: return c.c;
        case Gpu_Array::B: return c.B;
        case Gpu_Array::x: return c.x;
        case Gpu_Array::X: return c.X;
        case Gpu_Array::Y: return c.Y;
        }
      throw std::invalid_argument("Gpu_Solver: array");
    }

    Factor_Options xy_options() const
    {
      Factor_Options fo;
      if(!opt.correct_xy)
        fo.panel_refine = 0;
      return fo;
    }
    Factor_Options schur_options() const
    {
      Factor_Options fo;
      if(!opt.correct_schur)
        fo.panel_refine = 0;
      fo.correction_digits = opt.schur_correction_digits;
      fo.correction_limit = opt.schur_correction_limit;
      if(iteration == opt.corrupt_leaf_iteration)
        fo.corrupt_leaf_limbs = opt.corrupt_leaf_limbs;
      return fo;
    }

    Residue_Result residues()
    {
      Residue_Result r;
      stage("objectives");
      {
        Mat t(L, 1, N, 1);
        dev_mul(b.b, y.b, t.b, 0);
        r.dual_sum = mpf_total(t.b, 0);
      }
      for(auto &c : chunks)
        {
          Mat t(L, c.info.nb, c.info.P, 1);
          dev_mul(c.c.b, c.x.b, t.b, 0);
          r.primal_terms.push_back(mpf_total(t.b, 0));
        }

      stage("Cholesky X, Y + inverses");
      for(auto &c : chunks)
        {
          const auto &I = c.info;
          std::vector<size_t> ix(2 * I.nb);
          for(int q = 0; q < I.nb; ++q)
            for(int par = 0; par < 2; ++par)
              ix[par * I.nb + q] = 2 * size_t(I.blocks[q]) + par;
          auto valid = [&](int z, int i) {
            return i < I.n && (i % I.h) < c.hj[z];
          };
          factor_block(c.X, c.Xi, &c.XiT, xy_options(), conds.X, ix, valid);
          factor_block(c.Y, c.Yi, &c.YiT, xy_options(), conds.Y, ix, valid);
        }

      stage("pairings");
      for(auto &c : chunks)
        {
          Mat T = prod(c.XiT, c.Vx); // L_X^-1 Vx
          c.AXI = prod(T, T);
          mpf_lower_to_upper(c.AXI.b, 0);
          T.free();
          Mat YV = prod(c.Y, c.Vx); // Y Vx
          c.AY = prod(c.Vx, YV);
          mpf_lower_to_upper(c.AY.b, 0);
        }

      stage("dual residues");
      Host_Float dual_error;
      for(auto &c : chunks)
        {
          Mat d(L, c.info.nb, c.info.P, 1);
          mpf_seg_diag(c.AY.b, d.b, c.info.K, c.info.dim, 0);
          Mat By = gemv_B_y(c, y);
          Mat t = add(c.c, By, -1);
          By.free();
          c.dres = add(t, d, -1);
          dual_error = max_of(dual_error, mpf_max_abs(c.dres.b, 0));
        }
      dual_error = global_max(dual_error);
      r.dual_error = dual_error;
      // The floor of the dual residual that the explicit inverses of S leave only
      // matters once the residual has come down near it
      refine_system_now = opt.refine_system;
      if(opt.refine_system < 0)
        {
          const bool zero = dual_error.m.empty() || dual_error.e == -(int64_t(1) << 40);
          const int64_t top = zero ? 0 : dual_error.e + 32 * L - 2;
          refine_system_now = (!zero && top < opt.refine_system_below_log2) ? 1 : 0;
        }

      stage("primal residues P");
      Host_Float err_P;
      for(auto &c : chunks)
        {
          Mat Ax = cmws(c, c.x);
          c.Pres = add(Ax, c.X, -1);
          mpf_mask_blocks(c.Pres.b, c.info.h, c.d_hj, 0);
          err_P = max_of(err_P, mpf_max_abs(c.Pres.b, 0));
        }
      r.primal_error_P = global_max(err_P);

      stage("primal residues p");
      pvec = root_copy(b);
      for(auto &c : chunks)
        {
          Mat t = gemv_Bt(c, c.x);
          mpf_batch_add(pvec.b, t.b, -1, pvec.b, 0);
        }
      sum_ranks(pvec);
      r.primal_error_p = mpf_max_abs(pvec.b, 0);

      stage("X Y");
      for(auto &c : chunks)
        {
          const auto &I = c.info;
          Mat XY = prod(c.X, c.Y);
          mpf_mask_blocks(XY.b, I.h, c.d_hj, 0);
          dev_neg(XY.b, XY.b, 0);
          c.mXY = std::move(XY);
          Mat D(L, 2 * I.nb, I.n, 1);
          mpf_diagonal(c.mXY.b, D.b, 0);
          r.minus_xy_trace.push_back(mpf_total(D.b, 0));
        }
      stage(nullptr);
      return r;
    }

    void factor()
    {
      stage("Schur blocks");
      std::vector<Mat> S(chunks.size());
      for(size_t i = 0; i < chunks.size(); ++i)
        {
          auto &c = chunks[i];
          const auto &I = c.info;
          S[i] = Mat(L, I.nb, I.P, I.P);
          mpf_schur(c.AXI.b, c.AY.b, S[i].b, I.K, I.dim, c.d_Kj, 0);
          c.AXI.free();
          c.AY.free();
        }

      stage("Cholesky S + inverses");
      for(size_t i = 0; i < chunks.size(); ++i)
        {
          auto &c = chunks[i];
          const auto &I = c.info;
          std::vector<size_t> ix(I.nb);
          for(int q = 0; q < I.nb; ++q)
            ix[q] = size_t(I.blocks[q]);
          auto valid = [&](int z, int p) {
            return p < I.P && (p % I.K) < c.Kj[z];
          };
          factor_block(S[i], c.Si, &c.SiT, schur_options(), conds.schur, ix, valid);
          if(refine_system_now > 0)
            c.S = std::move(S[i]);
          else
            S[i].free();
        }

      stage("P = L^-1 B");
      int64_t rows_total = 0;
      for(auto &c : chunks)
        rows_total += int64_t(c.info.nb) * c.info.P;
      if(!gram && !world)
        {
          Gram_Options qopt = gopt;
          qopt.timing = opt.product_timing;
          gram = std::make_unique<Mpf_Gram>(N, L, D, rows_total, qopt,
                                            cudaStream_t(0));
        }
      for(auto &c : chunks)
        {
          const auto &I = c.info;
          const int64_t R = int64_t(I.nb) * I.P;
          if(!c.Pf)
            c.Pf = Mat(L, 1, R, N);
          if(opt.structured_p)
            eng(I.nb).solve_flat(c.SiT.b, B_of(c), c.Pf.b.m, c.Pf.b.e, opt.pq_guard);
          else
            eng(I.nb).mm_flat(c.SiT.b, B_of(c), c.Pf.b.m, c.Pf.b.e);
          B_done(c);
          eng(I.nb).check();
          c.SiT.free();
          if(gram) gram->scan(P_of(c).view()); // only takes maxima: needs one P at a time
        }
      form_and_factor_Q();
      stage(nullptr);
    }

    // Q = sum over chunks of P^T P and its factor: the one place that knows Q's
    // layout, with solve_Q (a distributed Q would replace both)
    void form_and_factor_Q()
    {
      if(world) { distributed_Q_factor(); return; }
      stage("Q = P^T P");
      for(auto &c : chunks)
        gram->add(P_of(c).view());
      if(opt.lean_p)
        for(auto &c : chunks)
          c.Pf.free();
      const int64_t packed = int64_t(N) * (N + 1) / 2;
      Mat Qp(L, 1, packed, 1), Q(L, 1, N, N);
      gram->finish(Qp.b.m, Qp.b.e);
      mpf_unpack_symmetric(Qp.b, Q.b, 0);
      Qp.free();
      if(want_Q_digest()) record_Q_digest(full_Q_tile_digests(Q.b,opt.q_tile));

      stage("Cholesky Q + inverse");
      // the optimized correction validated on real Q
      Factor_Options fo{true, 1, 0, -1, 1024, true};
      fo.crt_guard = opt.pq_guard;
      std::vector<Block_Cond> qc(1);
      std::vector<size_t> ix{0};
      auto valid = [&](int, int i) { return i < N; };
      factor_block(Q, Xq, nullptr, fo, qc, ix, valid);
      conds.Q = qc[0];
      if(opt.refine_dy > 0)
        Qfull = std::move(Q);
    }

    // Q^-1 v, through the explicit inverse of Q's factor, with opt.refine_dy
    // steps of iterative refinement against Q itself
    Mat solve_Q(const Mat &v)
    {
      if(world)
        {
          Mat out(L,1,N,1);
          distributed_q->solve(v.b,out.b,opt.refine_dy);
          record_distributed_Q_peak();
          return out;
        }
      auto through_inverse = [&](const Mat &r) {
        Mat t = gemv_square(Xq, r, false);
        return gemv_square(Xq, t, true);
      };
      Mat d = through_inverse(v);
      for(int k = 0; k < opt.refine_dy; ++k)
        {
          Mat Qd = gemv_square(Qfull, d, false);
          Mat r = add(v, Qd, -1);
          Qd.free();
          Mat e = through_inverse(r);
          mpf_batch_add(d.b, e.b, 1, d.b, 0);
        }
      return d;
    }

    Host_Float R_error(const Host_Float &mu)
    {
      stage("R error");
      Mat s = scalar(mu);
      Host_Float err;
      for(auto &c : chunks)
        {
          Mat R = make_copy(c.mXY);
          mpf_add_diag(R.b, s.b, c.info.h, c.d_hj, 0);
          err = max_of(err, mpf_max_abs(R.b, 0));
        }
      stage(nullptr);
      return global_max(err);
    }

    void direction(const Host_Float &beta_mu, bool corrector)
    {
      const char *const s_R = corrector ? "corrector: R" : "predictor: R";
      const char *const s_Z = corrector ? "corrector: Z" : "predictor: Z";
      const char *const s_rhs
        = corrector ? "corrector: Schur rhs" : "predictor: Schur rhs";
      const char *const s_w
        = corrector ? "corrector: w, dy rhs" : "predictor: w, dy rhs";
      const char *const s_dy = corrector ? "corrector: dy" : "predictor: dy";
      const char *const s_dx = corrector ? "corrector: dx" : "predictor: dx";
      const char *const s_dX = corrector ? "corrector: dX" : "predictor: dX";
      const char *const s_dY = corrector ? "corrector: dY" : "predictor: dY";

      Mat bm = scalar(beta_mu);
      Mat rhs_y = root_copy(pvec);
      for(auto &c : chunks)
        {
          const auto &I = c.info;
          stage(s_R);
          Mat R;
          if(corrector)
            {
              Mat dXdY = prod(c.dX, c.dY);
              R = add(c.mXY, dXdY, -1);
            }
          else
            R = make_copy(c.mXY);
          mpf_add_diag(R.b, bm.b, I.h, c.d_hj, 0);

          stage(s_Z);
          Mat Z;
          {
            Mat PY = prod(c.Pres, c.Y); // Pres Y (Pres symmetric)
            Mat W = add(PY, R, -1);
            PY.free();
            Mat U = solve_X(c, W);
            W.free();
            Z = symmetrize(U);
          }

          stage(s_rhs);
          Mat rhs;
          {
            Mat ZV = prod(Z, c.Vx);
            Mat VZV = prod(c.Vx, ZV);
            Mat d(L, I.nb, I.P, 1);
            mpf_seg_diag(VZV.b, d.b, I.K, I.dim, 0);
            Mat zero = make_zero(L, I.nb, I.P, 1);
            Mat t = add(zero, c.dres, -1);
            rhs = add(t, d, -1);
          }

          stage(s_w);
          c.w = gemv_square(c.Si, rhs, false);
          if(refine_system_now > 0 && corrector)
            c.rhs = make_copy(rhs);
          Mat dyc = gemv_Pt(c, c.w);
          mpf_batch_add(rhs_y.b, dyc.b, -1, rhs_y.b, 0);
          c.R = std::move(R);
        }

      stage(s_dy);
      sum_ranks(rhs_y);
      dy = solve_Q(rhs_y);

      for(auto &c : chunks)
        {
          stage(s_dx);
          Mat Pdy = gemv_P(c, dy);
          Mat t = add(c.w, Pdy, 1);
          Pdy.free();
          c.dx = gemv_square(c.Si, t, true);
          t.free();
        }
      if(corrector)
        for(int k = 0; k < refine_system_now; ++k)
          refine_system_step();

      for(auto &c : chunks)
        {
          const auto &I = c.info;
          stage(s_dX);
          Mat Ad = cmws(c, c.dx);
          c.dX = add(Ad, c.Pres, 1);
          Ad.free();

          stage(s_dY);
          Mat dXY = prod(c.dX, c.Y); // dX Y (dX symmetric)
          Mat W = add(dXY, c.R, -1);
          dXY.free();
          Mat U = solve_X(c, W);
          W.free();
          Mat Sy = symmetrize(U);
          U.free();
          dev_neg(Sy.b, Sy.b, 0);
          mpf_mask_blocks(Sy.b, I.h, c.d_hj, 0);
          c.dY = std::move(Sy);
          c.R.free();
          c.w.free();
        }
      stage(nullptr);
    }

    // X^-1 W through the explicit inverse of X's factor, with opt.refine_x steps
    // of iterative refinement against X itself
    Mat solve_X(Chunk &c, const Mat &W)
    {
      auto through_inverse = [&](const Mat &V) {
        Mat T = prod(c.XiT, V); // L^-1 V
        return prod(c.Xi, T);   // L^-T L^-1 V
      };
      Mat U = through_inverse(W);
      for(int k = 0; k < opt.refine_x; ++k)
        {
          Mat XU = prod(c.X, U); // X symmetric; zero padded rows of W stay zero
          mpf_mask_blocks(XU.b, c.info.h, c.d_hj, 0);
          Mat r = add(W, XU, -1);
          XU.free();
          Mat e = through_inverse(r);
          r.free();
          mpf_batch_add(U.b, e.b, 1, U.b, 0);
        }
      return U;
    }

    // One step of iterative refinement on the whole Schur system of the corrector,
    //   S dx - B dy = rx,  B^T dx = ry  (SDPB's solve_schur_complement_equation).
    // r1 = rx - S dx + B dy, r2 = ry - B^T dx, and the same block solve on (r1, r2).
    void refine_system_step()
    {
      stage("corrector: refine the Schur system");
      Mat rhs2 = root_copy(pvec);
      for(auto &c : chunks)
        {
          Mat Sdx = gemv_square(c.S, c.dx, false);
          Mat Bdy = gemv_B_y(c, dy);
          Mat t = add(c.rhs, Sdx, -1);
          Sdx.free();
          Mat r1 = add(t, Bdy, 1);
          Bdy.free();
          t.free();
          Mat Btdx = gemv_Bt(c, c.dx);
          mpf_batch_add(rhs2.b, Btdx.b, -1, rhs2.b, 0);
          Btdx.free();
          c.w2 = gemv_square(c.Si, r1, false);
          Mat ptw = gemv_Pt(c, c.w2);
          mpf_batch_add(rhs2.b, ptw.b, -1, rhs2.b, 0);
        }
      sum_ranks(rhs2);
      Mat ddy = solve_Q(rhs2);
      for(auto &c : chunks)
        {
          Mat Pd = gemv_P(c, ddy);
          Mat t = add(c.w2, Pd, 1);
          Pd.free();
          c.w2.free();
          Mat dd = gemv_square(c.Si, t, true);
          mpf_batch_add(c.dx.b, dd.b, 1, c.dx.b, 0);
        }
      mpf_batch_add(dy.b, ddy.b, 1, dy.b, 0);
    }

    std::vector<Host_Float> frobenius_terms()
    {
      stage("Frobenius product");
      std::vector<Host_Float> out;
      for(auto &c : chunks)
        {
          const auto &I = c.info;
          Mat A = add(c.X, c.dX, 1), B2 = add(c.Y, c.dY, 1);
          mpf_mask_blocks(A.b, I.h, c.d_hj, 0);
          mpf_mask_blocks(B2.b, I.h, c.d_hj, 0);
          Mat t(L, A.b.batch, A.b.rows, A.b.cols);
          dev_mul(A.b, B2.b, t.b, 0);
          out.push_back(mpf_total(t.b, 0));
        }
      stage(nullptr);
      return out;
    }

    Min_Eig min_eig(bool dual, double gamma, const Refine_Fn &refine)
    {
      Min_Eig out;
      ld best = 0;
      const double EPS = std::ldexp(1.0, -52);
      auto step_from = [&](ld lam) -> ld {
        return lam > -ld(gamma) ? ld(1) : -ld(gamma) / lam;
      };
      for(auto &c : chunks)
        {
          const auto &I = c.info;
          stage(dual ? "step: L^-1 dY L^-T" : "step: L^-1 dX L^-T");
          const Mat &LiT = dual ? c.YiT : c.XiT;
          const Mat &D = dual ? c.dY : c.dX;
          Mat T = prod(LiT, D); // L^-1 dM
          Mat Tt = transpose(T);
          T.free();
          Mat M = prod(LiT, Tt); // L^-1 dM L^-T
          Tt.free();

          stage("step: to double, transfer");
          const int n = I.n, batch = M.b.batch;
          double *d_vals;
          GS_CHECK(cudaMallocAsync(&d_vals, size_t(batch) * n * n * 8, 0));
          std::vector<int64_t> scale(batch);
          mpf_to_double(M.b, d_vals, scale.data(), 0);
          std::vector<double> vals(size_t(batch) * n * n);
          GS_CHECK(cudaMemcpyAsync(vals.data(), d_vals, vals.size() * 8,
                                   cudaMemcpyDeviceToHost, 0));
          GS_CHECK(cudaStreamSynchronize(0));
          GS_CHECK(cudaFreeAsync(d_vals, 0));

          stage("step: host eigenvalues");
          std::vector<double> lam(batch, 0), bound(batch, 0);
          std::vector<int> ok(batch, 0);
#pragma omp parallel
          {
            std::vector<double> A(size_t(n) * n), w(n);
            std::vector<double> work(1 + 2 * n + 16);
            std::vector<int> iwork(1);
            int lwork = int(work.size()), liwork = 1;
#pragma omp for schedule(dynamic)
            for(int z = 0; z < batch; ++z)
              {
                if(scale[z] == 0)
                  {
                    bool all_zero = true;
                    const double *p = vals.data() + size_t(z) * n * n;
                    for(int64_t q = 0; q < int64_t(n) * n; ++q)
                      if(p[q] != 0)
                        {
                          all_zero = false;
                          break;
                        }
                    if(all_zero)
                      continue;
                  }
                ok[z] = 1;
                const double *p = vals.data() + size_t(z) * n * n;
                double fro = 0;
                for(int j = 0; j < n; ++j)
                  for(int i = 0; i < n; ++i)
                    {
                      const double v = 0.5 * (p[j * n + i] + p[i * n + j]);
                      A[size_t(j) * n + i] = v;
                      fro += v * v;
                    }
                int info = 0;
                const char jobz = 'N', uplo = 'L';
                dsyevd_(&jobz, &uplo, &n, A.data(), &n, w.data(), work.data(),
                        &lwork, iwork.data(), &liwork, &info);
                if(info != 0)
                  throw std::runtime_error("Gpu_Solver: dsyevd failed");
                lam[z] = w[0];
                bound[z] = (n + 2) * EPS * std::sqrt(fro);
              }
          }

          stage("step: refine");
          std::vector<int> flagged;
          for(int z = 0; z < batch; ++z)
            {
              if(!ok[z])
                continue;
              const ld l = std::ldexp(ld(lam[z]), int(scale[z]));
              const ld e = std::ldexp(ld(bound[z]), int(scale[z]));
              if(step_from(l - e) >= (1 - ld(opt.eig_tol)) * step_from(l + e))
                {
                  if(!out.any || l < best)
                    {
                      best = l;
                      out.any = true;
                      // the double as an integer mantissa times 2^(exp + scale)
                      int ex;
                      const double f = std::frexp(lam[z], &ex);
                      const long long mi = (long long)std::ldexp(f, 53);
                      out.value.m = {uint32_t(mi), uint32_t(mi >> 32)};
                      out.value.e = int64_t(ex) - 53 + scale[z];
                    }
                }
              else
                flagged.push_back(z);
            }
          if(!flagged.empty())
            {
              std::vector<int32_t> hm(size_t(L) * M.b.entries());
              std::vector<int64_t> he(M.b.entries());
              mpf_download(M.b, hm.data(), he.data(), 0);
              const int64_t per = int64_t(n) * n, tot = M.b.entries();
              for(int z : flagged)
                {
                  std::vector<int32_t> bm(size_t(L) * per);
                  std::vector<int64_t> be(per);
                  for(int64_t q = 0; q < per; ++q)
                    {
                      be[q] = he[z * per + q];
                      for(int t = 0; t < L; ++t)
                        bm[t * per + q] = hm[t * tot + z * per + q];
                    }
                  const Host_Float v = refine(n, bm, be);
                  ++out.refined;
                  // the long double value decides the minimum only
                  ld l = 0;
                  if(!v.m.empty())
                    {
                      const int k = int(v.m.size());
                      const bool neg = int32_t(v.m[k - 1]) < 0;
                      ld acc = 0;
                      for(int t = k - 1; t >= std::max(0, k - 4); --t)
                        {
                          const uint32_t w
                            = neg ? ~v.m[t] : v.m[t];
                          acc = std::ldexp(acc, 32) + ld(w);
                        }
                      acc = std::ldexp(acc, int(std::max(0, k - 4) * 32));
                      l = std::ldexp(neg ? -acc : acc, int(v.e));
                    }
                  if(!out.any || l < best)
                    {
                      best = l;
                      out.any = true;
                      out.value = v;
                    }
                }
            }
        }
      if(opt.perturb_step != 0 && iteration == opt.perturb_iteration && out.any
         && !dual)
        {
          // negative control of the comparison with SDPB
          const long long scale = 1LL << 53;
          const long long fac
            = (long long)std::llround((1.0 + opt.perturb_step) * double(scale));
          Host_Float v = out.value;
          if(v.m.size() == 2)
            {
              const __int128 mant
                = (__int128)(long long)(((unsigned long long)v.m[1] << 32) | v.m[0])
                  * fac;
              const __int128 shifted = mant >> 30;
              out.value.m = {uint32_t(shifted), uint32_t(shifted >> 32),
                             uint32_t(shifted >> 64), uint32_t(shifted >> 96)};
              out.value.e = v.e + 30 - 53;
              fprintf(stderr, "Gpu_Solver: step eigenvalue perturbed at iteration %d\n",
                      iteration);
            }
        }
      stage(nullptr);
      return out;
    }

    void update(const Host_Float &aP, const Host_Float &aD)
    {
      stage("update");
      Mat sP = scalar(aP), sD = scalar(aD);
      for(auto &c : chunks)
        {
          {
            Mat t(L, c.dx.b.batch, c.dx.b.rows, 1);
            dev_mul(c.dx.b, sP.b, t.b, 0);
            mpf_batch_add(c.x.b, t.b, 1, c.x.b, 0);
          }
          {
            Mat t(L, c.dX.b.batch, c.dX.b.rows, c.dX.b.cols);
            dev_mul(c.dX.b, sP.b, t.b, 0);
            mpf_batch_add(c.X.b, t.b, 1, c.X.b, 0);
          }
          {
            Mat t(L, c.dY.b.batch, c.dY.b.rows, c.dY.b.cols);
            dev_mul(c.dY.b, sD.b, t.b, 0);
            mpf_batch_add(c.Y.b, t.b, 1, c.Y.b, 0);
          }
          c.dx.free();
          c.dX.free();
          c.dY.free();
          c.Xi.free();
          c.XiT.free();
          c.Yi.free();
          c.YiT.free();
          c.Si.free();
          c.Pf.free();
          c.S.free();
          c.rhs.free();
          c.dres.free();
          c.Pres.free();
          c.mXY.free();
        }
      Mat t(L, 1, N, 1);
      dev_mul(dy.b, sD.b, t.b, 0);
      mpf_batch_add(y.b, t.b, 1, y.b, 0);
      dy.free();
      Xq.free();
      Qfull.free();
      if(distributed_q)
        {
          distributed_profiles.push_back(distributed_q->profile_json());
          peak_bytes=std::max(peak_bytes,double(distributed_q->peak_device_bytes()));
          distributed_q.reset();
        }
      stage(nullptr);
    }

    void end_iteration()
    {
      stage("CRT checks");
      for(auto &kv : engines)
        kv.second->check();
      stage(nullptr);
      last_iter_s = iter_s;
      iter_s.clear();
      last_iter_payload=iter_payload;
      iter_payload.clear();
      last_iter_comm_s=iter_comm_s;
      iter_comm_s.clear();
      ++iteration;
    }

    std::string profile_json() const
    {
      std::ostringstream os;
      os.precision(6);
      auto put = [&](const auto &m) {
        os << "{";
        bool first = true;
        for(const auto &kv : m)
          {
            os << (first ? "" : ", ") << "\"" << kv.first << "\": " << kv.second;
            first = false;
          }
        os << "}";
      };
      std::map<std::string, double> product_ms;
      if(opt.product_timing)
        {
          for(const auto &e : engines)
            e.second->collect_product_timings(product_ms);
          if(gram)
            {
              for(const auto &kv : gram->timings())
                product_ms["Q: " + kv.first] += kv.second;
              for(const auto &kv : gram->gram().timings())
                product_ms["Q: " + kv.first] += kv.second;
            }
        }
      size_t evictions = 0, product_bytes = workspace ? workspace->bytes() : 0;
      for(const auto &e : engines)
        evictions += e.second->cache_evictions();
      os << "{\"iterations\": " << iteration << ", \"peak_device_bytes\": "
         << peak_bytes << ", \"cache_evictions\": " << evictions
         << ", \"workspace_bytes\": " << product_bytes << ", \"total_s\": ";
      put(total_s);
      os << ", \"last_iteration_s\": ";
      put(last_iter_s);
      os << ", \"device_high_water_gb_by_stage\": ";
      put(stage_used_gb);
      os << ", \"host_peak_rss_gb_by_stage\": ";
      put(stage_host_gb);
      os << ", \"product_phases_ms\": ";
      put(product_ms);
      if(world)
        {
          os << ", \"collective_received_payload_by_stage\": ";
          put(total_payload);
          os << ", \"last_iteration_collective_received_payload_by_stage\": ";
          put(last_iter_payload);
          os << ", \"rank\": " << world->rank()
             << ", \"G\": " << world->size()
             << ", \"collective_received_payload\": " << world->received_bytes();
          if(world->timing())
            {
              os << ", \"collective_seconds_by_stage\": ";
              put(total_comm_s);
              os << ", \"last_iteration_collective_seconds_by_stage\": ";
              put(last_iter_comm_s);
              os << ", \"collective_seconds\": {\"broadcast\": " << world->broadcast_seconds()
                 << ", \"max_u64\": " << world->max_seconds()
                 << ", \"reduce_scatter\": " << world->reduce_scatter_seconds()
                 << ", \"allgather\": " << world->allgather_seconds() << "}";
            }
          os << ", \"distributed_q_iterations\": [";
          for(size_t i=0;i<distributed_profiles.size();++i)
            os << (i?",":"") << distributed_profiles[i];
          os << "]";
        }
      os << "}\n";
      return os.str();
    }
  };

  Gpu_Solver::Gpu_Solver(const Gpu_Solver_Options &options, int N,
                         const std::vector<Block_Shape> &blocks)
      : p_(std::make_unique<Impl>())
  {
    p_->opt = options;
    p_->N = N;
    p_->L = options.L;
    p_->D = options.digits;
    p_->shapes = blocks;
    p_->gopt.prime_group = 32;
    p_->gopt.fused_residues = true;
    p_->gopt.product_timing = options.product_timing;
    p_->gopt.crt_guard_bits = options.crt_guard;
    if(options.world.size>1)
      p_->world=std::make_unique<Gpu_Collective>(options.world);
    p_->workspace = std::make_unique<Rns_Workspace>(cudaStream_t(0));
    p_->gopt.workspace = p_->workspace.get();
    g_workspace = p_->workspace.get();
    p_->build_chunks();
    p_->b = Mat(p_->L, 1, N, 1);
    p_->y = Mat(p_->L, 1, N, 1);
  }

  Gpu_Solver::~Gpu_Solver()
  {
    if(p_)
      for(auto &c : p_->chunks)
        {
          cudaFree(c.d_hj);
          cudaFree(c.d_Kj);
        }
  }

  int Gpu_Solver::num_chunks() const { return int(p_->chunks.size()); }
  const Chunk_Info &Gpu_Solver::chunk(int i) const { return p_->chunks.at(i).info; }
  int64_t Gpu_Solver::entries(int i, Gpu_Array a) const
  {
    return p_->array(i, a).b.entries();
  }
  void Gpu_Solver::upload(int i, Gpu_Array a, const int32_t *m, const int64_t *e)
  {
    if(a == Gpu_Array::B && p_->opt.sparse_b)
      p_->upload_sparse_B(p_->chunks.at(i), m, e);
    else
      mpf_upload(m, e, p_->array(i, a).b, 0);
  }
  void Gpu_Solver::download(int i, Gpu_Array a, int32_t *m, int64_t *e)
  {
    mpf_download(p_->array(i, a).b, m, e, 0);
  }
  void Gpu_Solver::upload_b(const int32_t *m, const int64_t *e)
  {
    mpf_upload(m, e, p_->b.b, 0);
  }
  void Gpu_Solver::upload_y(const int32_t *m, const int64_t *e)
  {
    mpf_upload(m, e, p_->y.b, 0);
  }
  void Gpu_Solver::download_y(int32_t *m, int64_t *e)
  {
    mpf_download(p_->y.b, m, e, 0);
  }
  void Gpu_Solver::finish_load()
  {
    auto &P = *p_;
    const int L = P.L;
    for(auto &c : P.chunks)
      {
        const auto &I = c.info;
        c.VT = P.transpose(c.V);
        c.VTrep = Mat(L, 2 * I.nb * c.nseg, I.K, I.h);
        mpf_repeat(c.VT.b, c.VTrep.b, 0);
        c.Vx = Mat(L, 2 * I.nb, I.n, I.dim * I.K);
        mpf_kron_identity(c.V.b, c.Vx.b, I.dim, 0);
        c.V.free();
      }
    GS_CHECK(cudaStreamSynchronize(0));
  }
  void Gpu_Solver::end_run()
  {
    auto &P = *p_;
    for(auto &c : P.chunks)
      {
        c.Xi.free();
        c.XiT.free();
        c.Yi.free();
        c.YiT.free();
        c.AXI.free();
        c.AY.free();
        c.dres.free();
        c.Pres.free();
        c.mXY.free();
      }
    P.pvec.free();
    GS_CHECK(cudaStreamSynchronize(0));
  }
  Residue_Result Gpu_Solver::residues() { return p_->residues(); }
  void Gpu_Solver::factor() { p_->factor(); }
  Host_Float Gpu_Solver::R_error(const Host_Float &mu)
  {
    return p_->R_error(mu);
  }
  void Gpu_Solver::direction(const Host_Float &beta_mu, bool corrector)
  {
    p_->direction(beta_mu, corrector);
  }
  std::vector<Host_Float> Gpu_Solver::frobenius_terms()
  {
    return p_->frobenius_terms();
  }
  Min_Eig Gpu_Solver::min_eig(bool dual, double gamma, const Refine_Fn &refine)
  {
    return p_->min_eig(dual, gamma, refine);
  }
  void Gpu_Solver::update(const Host_Float &aP, const Host_Float &aD)
  {
    p_->update(aP, aD);
  }
  const Cond_Table &Gpu_Solver::cond() const { return p_->conds; }
  void Gpu_Solver::end_iteration() { p_->end_iteration(); }
  std::string Gpu_Solver::profile_json() const { return p_->profile_json(); }
  double Gpu_Solver::peak_device_bytes() const { return p_->peak_bytes; }
}
