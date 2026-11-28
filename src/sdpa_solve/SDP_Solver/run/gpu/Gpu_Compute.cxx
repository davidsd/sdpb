#include "Gpu_Compute.hxx"
#include "gpu/solver/Device_Convert.hxx"
#include "gpu/rns/Limb_Sizes.hxx"
#include "gpu/solver/Gpu_Solver.hxx"
#include "sdp_solve/SDP_Solver/run/step/centering_parameters.hxx"
#include "sdpa_solve/SDP_Solver.hxx"
#include "sdpa_solve/parse_sdpa/parse_sdpa.hxx"
#include "sdpb_util/assert.hxx"

#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <tuple>
#include <numeric>

using namespace sdpb_gpu;

namespace Sdpb::Sdpa
{
  namespace
  {
    // Lower-triangle entries of a chunk read by one slot of the reader, in
    // file order. Index -1 denotes F_0.
    struct F_Piece
    {
      std::vector<int32_t> p;
      std::vector<uint32_t> pos;
      std::vector<uint64_t> order; // place in the file
      std::vector<int32_t> m;      // L limbs per entry, entry by entry
      std::vector<int64_t> e;
    };
  }

  struct Gpu_Compute::Impl
  {
    std::unique_ptr<Gpu_Solver> gs;
    int L, m;
    std::vector<int> dims, local_of, chunk_of, q_of;
    Residue_Result last;
    // Seconds spent reading and converting F, and packing and uploading it
    double read_F_s = 0, upload_F_s = 0;

    template <class Fn> void fill(int chunk, Gpu_Array a, const Fn &fn)
    {
      const int64_t tot = gs->entries(chunk, a);
      std::vector<int32_t> mm(size_t(L) * tot);
      std::vector<int64_t> e(tot);
      parallel_for(tot, [&](int64_t q) { fn(q, mm.data(), tot, e.data()); });
      gs->upload(chunk, a, mm.data(), e.data());
    }

    void upload_vector(const Primal_Dist_Vector &v, bool objective)
    {
      std::vector<int32_t> mm(size_t(L) * m);
      std::vector<int64_t> e(m);
      const auto &M = v.LockedMatrix();
      for(int i = 0; i < m; ++i)
        to_device(M.CRef(i, 0), L, mm.data(), m, i, e.data());
      if(objective)
        gs->upload_b(mm.data(), e.data());
      else
        gs->upload_y(mm.data(), e.data());
    }

    void upload_X_Y(int ci, const SDP_Solver &solver)
    {
      const Chunk_Info &I = gs->chunk(ci);
      const int n = I.n;
      auto fill_matrix = [&](Gpu_Array which, const Block_Diagonal_Matrix &M) {
        fill(ci, which, [&](int64_t q, int32_t *mm, int64_t tot, int64_t *e) {
          const int64_t per = int64_t(n) * n;
          const int z = int(q / per);
          const int64_t rem = q % per;
          const int i = int(rem % n), j = int(rem / n);
          const int b = I.blocks[z];
          if(i < dims[b] && j < dims[b])
            to_device(M.blocks[local_of[b]].LockedMatrix().CRef(i, j), L, mm, tot,
                      q, e);
          else if(i == j)
            one_device(L, mm, tot, q, e);
          else
            zero_device(L, mm, tot, q, e);
        });
      };
      fill_matrix(Gpu_Array::X, solver.X);
      fill_matrix(Gpu_Array::Y, solver.Y);
    }

    // Sparse inputs may repeat an entry. Match the CPU reader's last-write
    // rule, including an explicit zero that cancels an earlier value.
    void load_F(const std::filesystem::path &path)
    {
      const auto t_read = std::chrono::steady_clock::now();
      const int nc = gs->num_chunks();
      const size_t slots = sdpa_visit_slots();
      // pieces[slot][chunk]
      std::vector<std::vector<F_Piece>> pieces(slots, std::vector<F_Piece>(nc));
      std::vector<std::vector<int32_t>> limbs(slots, std::vector<int32_t>(L));
      visit_sdpa_lower(path, [&](size_t slot, uint64_t order, size_t block,
                                 size_t matrix, size_t i, size_t j,
                                 const El::BigFloat &value) {
        const int ci = chunk_of.at(block);
        if(ci < 0)
          return;
        const int64_t n = gs->chunk(ci).n;
        const auto pos = int64_t(q_of[block]) * n * n + int64_t(j) * n + i;
        ASSERT(pos <= std::numeric_limits<uint32_t>::max(), "GPU chunk is too large");
        auto &part = pieces.at(slot)[ci];
        auto &l = limbs[slot];
        int64_t exponent;
        to_device(value, L, l.data(), 1, 0, &exponent);
        part.p.push_back(int32_t(matrix) - 1);
        part.pos.push_back(uint32_t(pos));
        part.order.push_back(order);
        part.m.insert(part.m.end(), l.begin(), l.end());
        part.e.push_back(exponent);
      });
      const auto t_upload = std::chrono::steady_clock::now();
      read_F_s = std::chrono::duration<double>(t_upload - t_read).count();
      for(int ci = 0; ci < nc; ++ci)
        {
          const Chunk_Info &I = gs->chunk(ci);
          const int64_t n = I.n;
          // Entries as (slot, index), sorted by matrix, position and file order
          using Ref = std::pair<uint32_t, uint32_t>;
          std::vector<Ref> refs;
          for(size_t s = 0; s < slots; ++s)
            {
              const size_t count = pieces[s][ci].p.size();
              ASSERT(count <= std::numeric_limits<uint32_t>::max());
              for(size_t k = 0; k < count; ++k)
                refs.emplace_back(uint32_t(s), uint32_t(k));
            }
          auto key = [&](const Ref &r) {
            const F_Piece &P = pieces[r.first][ci];
            return std::make_tuple(P.p[r.second], P.pos[r.second],
                                   P.order[r.second]);
          };
          std::sort(refs.begin(), refs.end(),
                    [&](const Ref &a, const Ref &b) { return key(a) < key(b); });
          // The last line of each position; no explicit zeros
          auto same_place = [&](const Ref &a, const Ref &b) {
            return std::get<0>(key(a)) == std::get<0>(key(b))
                   && std::get<1>(key(a)) == std::get<1>(key(b));
          };
          std::vector<Ref> keep;
          for(size_t t = 0; t < refs.size(); ++t)
            {
              if(t + 1 < refs.size() && same_place(refs[t], refs[t + 1]))
                continue;
              if(pieces[refs[t].first][ci].e[refs[t].second] != MPF_ZERO_EXP)
                keep.push_back(refs[t]);
            }
          refs.clear();
          refs.shrink_to_fit();
          size_t nF0 = 0;
          while(nF0 < keep.size() && pieces[keep[nF0].first][ci].p[keep[nF0].second] < 0)
            ++nF0;
          {
            const int64_t tot = int64_t(I.nb) * n * n;
            std::vector<int32_t> mm(size_t(L) * tot, 0);
            std::vector<int64_t> e(tot, MPF_ZERO_EXP);
            for(size_t t = 0; t < nF0; ++t)
              {
                const F_Piece &part = pieces[keep[t].first][ci];
                const size_t k = keep[t].second;
                const int64_t a = part.pos[k], base = a / (n * n) * (n * n), rem = a % (n * n);
                const int64_t b = base + (rem % n) * n + rem / n;
                for(int w = 0; w < L; ++w)
                  {
                    mm[w * tot + a] = part.m[size_t(L) * k + w];
                    mm[w * tot + b] = part.m[size_t(L) * k + w];
                  }
                e[a] = part.e[k];
                e[b] = part.e[k];
              }
            gs->upload(ci, Gpu_Array::c, mm.data(), e.data());
          }
          const int64_t nnz = int64_t(keep.size() - nF0);
          std::vector<int32_t> pv(nnz), mm(size_t(L) * nnz);
          std::vector<uint32_t> posv(nnz);
          std::vector<int64_t> e(nnz);
          parallel_for(nnz, [&](int64_t t) {
            const F_Piece &part = pieces[keep[nF0 + t].first][ci];
            const size_t k = keep[nF0 + t].second;
            pv[t] = part.p[k];
            posv[t] = part.pos[k];
            e[t] = part.e[k];
            for(int w = 0; w < L; ++w)
              mm[w * nnz + t] = part.m[size_t(L) * k + w];
          });
          gs->sdpa_upload_F(ci, nnz, pv.data(), posv.data(), mm.data(), e.data());
          for(size_t s = 0; s < slots; ++s)
            pieces[s][ci] = F_Piece();
        }
      upload_F_s = std::chrono::duration<double>(std::chrono::steady_clock::now()
                                                 - t_upload)
                     .count();
    }
  };

  Gpu_Compute::Gpu_Compute(const Solver_Parameters &parameters,
                           const Block_Info &block_info, const SDP &sdp,
                           const SDP_Solver &solver)
      : p_(std::make_unique<Impl>())
  {
    ASSERT(El::mpi::Size(block_info.mpi_comm.value) == 1,
           "--gpu requires complete blocks on each MPI rank");
    auto &P = *p_;
    // limbs for --precision, rounded up to a size with compiled kernels
    P.L = sdpb_gpu::rounded_L(int((parameters.precision + 2 + 31) / 32) + 1);
    ASSERT(P.L > 0, "--gpu: --precision above the sizes in gpu/rns/Limb_Sizes.hxx");
    ASSERT(block_info.primal_dimension <= size_t(std::numeric_limits<int>::max()));
    P.m = int(block_info.primal_dimension);
    const int num_blocks = int(block_info.block_dimensions.size());
    P.dims.resize(num_blocks);
    P.local_of.assign(num_blocks, -1);
    P.chunk_of.assign(num_blocks, -1);
    P.q_of.assign(num_blocks, -1);
    for(int j = 0; j < num_blocks; ++j)
      P.dims[j] = int(block_info.block_dimensions[j]);
    for(size_t i = 0; i < block_info.block_indices.size(); ++i)
      P.local_of[block_info.block_indices[i]] = int(i);

    Gpu_Solver_Options o;
    o.world.rank = El::mpi::Rank();
    o.world.size = El::mpi::Size();
    o.world.device = -1;
    o.world.local_rank = block_info.node_comm.Rank();
    o.world.host_broadcast = [](void *p, size_t n, int root) {
      El::mpi::Broadcast(static_cast<unsigned char *>(p), int(n), root,
                         El::mpi::COMM_WORLD);
    };
    o.q_tile = int(env_long("SDPB_GPU_Q_TILE", 512));
    o.q_row_panel = int(env_long("SDPB_GPU_Q_ROW_PANEL", 8192));
    o.q_wide_residual = env_long("SDPB_GPU_Q_WIDE_RESIDUAL", 0) != 0;
    o.q_single_pass = env_long("SDPB_GPU_Q_SINGLE_PASS", 1) != 0;
    o.structured_q = env_long("SDPB_GPU_STRUCTURED_Q", 1) != 0;
    o.q_base = int(env_long("SDPB_GPU_Q_BASE", 2048));
    o.L = P.L;
    o.digits = int(env_long("SDPB_GPU_DIGITS", 4 * P.L));
    o.eig_tol = env_double("SDPB_GPU_EIG_TOL", 1e-10);
    o.correct_xy = env_long("SDPB_GPU_CORRECT_XY", 1) != 0;
    o.refine_dy = int(env_long("SDPB_GPU_REFINE_DY", 1));
    o.refine_x = int(env_long("SDPB_GPU_REFINE_X", 0));
    // Keep complete products for ill-conditioned SDPA systems.
    o.crt_guard = int(env_long("SDPB_GPU_CRT_GUARD", 0));
    o.pq_guard = int(env_long("SDPB_GPU_PQ_GUARD", 0));
    o.max_chunk_blocks = int(env_long("SDPB_GPU_CHUNK_BLOCKS", 0));
    o.max_chunk_rows = env_long("SDPB_GPU_CHUNK_ROWS", 0);
    // The streamed Schur assembly uses one-gigabyte bands.
    setenv("SDPB_GPU_PANEL", "4096", 0);
    o.product_timing = env_long("SDPB_GPU_PRODUCT_TIMING", 0) != 0;
    o.perturb_step = env_double("SDPB_GPU_PERTURB_STEP", 0);
    o.perturb_iteration = int(env_long("SDPB_GPU_PERTURB_ITERATION", 0));
    const auto t_init = std::chrono::steady_clock::now();
    P.gs = std::make_unique<Gpu_Solver>(o, P.m, P.dims);
    const double init_s = std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t_init)
                            .count();

    for(int ci = 0; ci < P.gs->num_chunks(); ++ci)
      {
        const Chunk_Info &I = P.gs->chunk(ci);
        for(int q = 0; q < I.nb; ++q)
          {
            ASSERT(P.local_of[I.blocks[q]] >= 0, "--gpu: chunk ", ci, " has block ",
                   I.blocks[q], ", which this rank doesn't hold");
            P.chunk_of[I.blocks[q]] = ci;
            P.q_of[I.blocks[q]] = q;
          }
      }
    const auto t0 = std::chrono::steady_clock::now();
    const auto since = [](const std::chrono::steady_clock::time_point &t) {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - t)
        .count();
    };
    P.upload_vector(sdp.primal_objective_c, true);
    P.upload_vector(solver.x, false);
    P.load_F(sdp.gpu_input);
    const auto t_xy = std::chrono::steady_clock::now();
    for(int ci = 0; ci < P.gs->num_chunks(); ++ci)
      P.upload_X_Y(ci, solver);
    const double xy_s = since(t_xy);
    const auto t_finish = std::chrono::steady_clock::now();
    P.gs->finish_load();
    const double finish_s = since(t_finish);
    if(El::mpi::Rank() == 0)
      {
        std::cout << "Gpu_Compute (SDPA): " << P.gs->num_chunks()
                  << " chunks, L = " << P.L << ", m = " << P.m << ", loaded in "
                  << since(t0) << " s (device setup " << init_s
                  << " s before; read F " << P.read_F_s << " s, upload F "
                  << P.upload_F_s << " s, X and Y " << xy_s << " s, finish "
                  << finish_s << " s)\n";
        for(int ci = 0; ci < P.gs->num_chunks(); ++ci)
          std::cout << "  chunk " << ci << ": " << P.gs->chunk(ci).nb
                    << " blocks, n " << P.gs->chunk(ci).n << "\n";
        std::cout << std::flush;
      }
  }

  Gpu_Compute::~Gpu_Compute() = default;

  void Gpu_Compute::residues(const SDP &sdp, SDP_Solver &solver)
  {
    auto &P = *p_;
    P.last = P.gs->sdpa_residues();
    El::BigFloat dsum(0);
    for(const auto &t : P.last.primal_terms)
      dsum += from_host(t);
    if(El::mpi::Size() > 1)
      dsum = El::mpi::AllReduce(dsum, El::mpi::SUM, El::mpi::COMM_WORLD);
    solver.primal_objective = sdp.objective_const + from_host(P.last.dual_sum);
    solver.dual_objective = sdp.objective_const + dsum;
    solver.duality_gap
      = El::Abs(solver.primal_objective - solver.dual_objective)
        / El::Max(El::Abs(solver.primal_objective)
                    + El::Abs(solver.dual_objective),
                  El::BigFloat(1));
    solver.primal_error = from_host(P.last.primal_error_P);
    solver.dual_error = from_host(P.last.dual_error);
  }

  void Gpu_Compute::step(const Solver_Parameters &parameters,
                         size_t total_psd_rows,
                         bool is_primal_and_dual_feasible, SDP_Solver &solver,
                         El::BigFloat &mu, El::BigFloat &beta_corrector,
                         El::BigFloat &primal_step_length,
                         El::BigFloat &dual_step_length, bool &terminate_now,
                         El::BigFloat &S_cond_number,
                         El::BigFloat &max_block_cond_number,
                         std::string &max_block_cond_number_name)
  {
    auto &P = *p_;
    const int L = P.L;
    El::BigFloat trace(0);
    for(const auto &t : P.last.minus_xy_trace)
      trace += from_host(t);
    if(El::mpi::Size() > 1)
      trace = El::mpi::AllReduce(trace, El::mpi::SUM, El::mpi::COMM_WORLD);
    mu = -trace / total_psd_rows;
    if(mu > parameters.max_complementarity)
      {
        terminate_now = true;
        return;
      }
    P.gs->sdpa_factor();
    solver.R_error = from_host(P.gs->R_error(to_host(mu, L)));

    const El::BigFloat beta_predictor
      = predictor_centering_parameter(parameters, is_primal_and_dual_feasible);
    P.gs->sdpa_direction(to_host(beta_predictor * mu, L), false);

    // corrector_centering_parameter.hxx
    El::BigFloat frobenius(0);
    for(const auto &t : P.gs->frobenius_terms())
      frobenius += from_host(t);
    if(El::mpi::Size() > 1)
      frobenius = El::mpi::AllReduce(frobenius, El::mpi::SUM, El::mpi::COMM_WORLD);
    const El::BigFloat r = frobenius / (mu * total_psd_rows);
    const El::BigFloat beta = r < 1 ? r * r : r;
    beta_corrector
      = is_primal_and_dual_feasible
          ? El::Min(El::Max(parameters.feasible_centering_parameter, beta),
                    El::BigFloat(1))
          : El::Max(parameters.infeasible_centering_parameter, beta);
    P.gs->sdpa_direction(to_host(beta_corrector * mu, L), true);

    // update_cond_numbers.hxx, same order of blocks
    const Cond_Table &C = P.gs->cond();
    auto cond_of = [&](const Block_Cond &c) {
      if(c.max.m.empty() || c.min.m.empty())
        return El::BigFloat(0);
      const El::BigFloat ratio = from_host(c.max) / from_host(c.min);
      return ratio * ratio;
    };
    S_cond_number = cond_of(C.Q);
    max_block_cond_number = 0;
    max_block_cond_number_name = "";
    int max_condition_order = std::numeric_limits<int>::max();
    for(size_t j = 0; j < C.X.size(); ++j)
      {
        const El::BigFloat x = cond_of(C.X[j]);
        if(max_block_cond_number < x)
          {
            max_block_cond_number = x;
            max_condition_order = int(2 * j);
            max_block_cond_number_name = El::BuildString("X_cholesky.block_", j);
          }
        const El::BigFloat y = cond_of(C.Y[j]);
        if(max_block_cond_number < y)
          {
            max_block_cond_number = y;
            max_condition_order = int(2 * j + 1);
            max_block_cond_number_name = El::BuildString("Y_cholesky.block_", j);
          }
      }
    if(El::mpi::Size() > 1)
      {
        const auto local = max_block_cond_number;
        max_block_cond_number
          = El::mpi::AllReduce(local, El::mpi::MAX, El::mpi::COMM_WORLD);
        // as the one-rank traversal on exact ties, whoever owns the block
        const int order = local == max_block_cond_number
                            ? max_condition_order
                            : std::numeric_limits<int>::max();
        const int first = El::mpi::AllReduce(order, El::mpi::MIN, El::mpi::COMM_WORLD);
        int root = order == first ? El::mpi::Rank() : El::mpi::Size();
        root = El::mpi::AllReduce(root, El::mpi::MIN, El::mpi::COMM_WORLD);
        int length = int(max_block_cond_number_name.size());
        El::mpi::Broadcast(&length, 1, root, El::mpi::COMM_WORLD);
        max_block_cond_number_name.resize(length);
        El::mpi::Broadcast(
          reinterpret_cast<unsigned char *>(max_block_cond_number_name.data()),
          length, root, El::mpi::COMM_WORLD);
      }

    // a block that the float64 eigenvalue cannot resolve: SDPB's multiprecision
    // eigenvalue (min_eigenvalue.cxx) on the host
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
    auto step_of = [&](const Min_Eig &me) {
      if(!me.any)
        return El::BigFloat(1);
      const El::BigFloat lam = from_host(me.value);
      if(lam > -parameters.step_length_reduction)
        return El::BigFloat(1);
      return -parameters.step_length_reduction / lam;
    };
    primal_step_length = step_of(lamP);
    dual_step_length = step_of(lamD);
    if(El::mpi::Size() > 1)
      {
        primal_step_length
          = El::mpi::AllReduce(primal_step_length, El::mpi::MIN, El::mpi::COMM_WORLD);
        dual_step_length
          = El::mpi::AllReduce(dual_step_length, El::mpi::MIN, El::mpi::COMM_WORLD);
      }
    if(lamP.refined + lamD.refined > 0)
      std::cerr << "Gpu_Compute: " << lamP.refined + lamD.refined
                << " blocks' step eigenvalues refined in multiprecision\n";
    if(is_primal_and_dual_feasible)
      {
        primal_step_length = El::Min(primal_step_length, dual_step_length);
        dual_step_length = primal_step_length;
      }
    P.gs->sdpa_update(to_host(primal_step_length, L), to_host(dual_step_length, L));
    P.gs->end_iteration();
  }

  void Gpu_Compute::get_iterate(SDP_Solver &solver)
  {
    auto &P = *p_;
    const int L = P.L;
    auto bf = [&](const std::vector<int32_t> &m, const std::vector<int64_t> &e,
                  int64_t tot, int64_t q) {
      return e[q] == MPF_ZERO_EXP
               ? El::BigFloat(0)
               : from_limbs(reinterpret_cast<const uint32_t *>(m.data()) + q, L,
                            tot, e[q]);
    };
    {
      std::vector<int32_t> m(size_t(L) * P.m);
      std::vector<int64_t> e(P.m);
      P.gs->download_y(m.data(), e.data());
      auto &x = solver.x.Matrix();
      for(int i = 0; i < P.m; ++i)
        x.Set(i, 0, bf(m, e, P.m, i));
    }
    std::vector<int32_t> m;
    std::vector<int64_t> e;
    for(int ci = 0; ci < P.gs->num_chunks(); ++ci)
      {
        const Chunk_Info &I = P.gs->chunk(ci);
        const int n = I.n;
        auto get_matrix = [&](Gpu_Array which, Block_Diagonal_Matrix &M) {
          const int64_t tot = P.gs->entries(ci, which);
          m.resize(size_t(L) * tot);
          e.resize(tot);
          P.gs->download(ci, which, m.data(), e.data());
          const int64_t per = int64_t(n) * n;
          for(int z = 0; z < I.nb; ++z)
            {
              const int b = I.blocks[z];
              auto &Mb = M.blocks[P.local_of[b]].Matrix();
              parallel_for(per, [&](int64_t q) {
                const int i = int(q % n), j = int(q / n);
                if(i < P.dims[b] && j < P.dims[b])
                  Mb.Set(i, j, bf(m, e, tot, z * per + q));
              });
            }
        };
        get_matrix(Gpu_Array::X, solver.X);
        get_matrix(Gpu_Array::Y, solver.Y);
      }
  }

  void Gpu_Compute::write_profile(const std::filesystem::path &path) const
  {
    const auto dest
      = El::mpi::Rank() == 0
          ? path
          : std::filesystem::path(path.string() + ".rank"
                                  + std::to_string(El::mpi::Rank()));
    std::ofstream os(dest);
    os << p_->gs->profile_json();
  }
}
