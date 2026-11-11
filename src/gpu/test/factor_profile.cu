// factor + inverse + X B on real dumps
// (Schur blocks or Q, bench/sdpvibe_dump_schur.py's layout), by stage.
//   factor_profile --dir <dump> --chunk 0 --nb 36 --P 945 --N 3406 --reps 6
//   --profile-reps 3 [Factor_Options flags]
// Three separate runs of the same work, each with a fresh Block_Factor:
//   TIME: nothing synchronized but the ends of the factor and of the solve,
//   STAGE: a synchronization after every Block_Factor stage,
//   KERNEL: also after every stage inside every product.
// The first repetition of each builds the product cache and is not counted.
#include "../linalg/Block_Factor.hxx"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

using namespace sdpb_gpu;

namespace
{
  const char *arg(int argc, char **argv, const char *name, const char *dflt)
  {
    for(int i = 1; i + 1 < argc; ++i)
      if(!strcmp(argv[i], name))
        return argv[i + 1];
    return dflt;
  }
  template <class T> std::vector<T> slurp(const std::string &path, size_t count)
  {
    std::vector<T> v(count);
    FILE *f = fopen(path.c_str(), "rb");
    if(!f || fread(v.data(), sizeof(T), count, f) != count)
      throw std::runtime_error("read " + path);
    fclose(f);
    return v;
  }
  double seconds(std::chrono::steady_clock::time_point a,
                 std::chrono::steady_clock::time_point b)
  {
    return std::chrono::duration<double>(b - a).count();
  }
  // the stream-ordered pool's used high-water mark since the last call, bytes
  double pool_high()
  {
    cudaMemPool_t pool;
    cudaDeviceGetDefaultMemPool(&pool, 0);
    size_t high = 0, zero = 0;
    cudaMemPoolGetAttribute(pool, cudaMemPoolAttrUsedMemHigh, &high);
    cudaMemPoolSetAttribute(pool, cudaMemPoolAttrUsedMemHigh, &zero);
    return double(high);
  }
  // the same for what the pool holds from the device, used or not
  double pool_reserved_high()
  {
    cudaMemPool_t pool;
    cudaDeviceGetDefaultMemPool(&pool, 0);
    size_t high = 0, zero = 0;
    cudaMemPoolGetAttribute(pool, cudaMemPoolAttrReservedMemHigh, &high);
    cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReservedMemHigh, &zero);
    return double(high);
  }
  double host_high()
  {
    struct rusage u;
    getrusage(RUSAGE_SELF, &u);
    return double(u.ru_maxrss) * 1024;
  }
  double median(std::vector<double> v)
  {
    if(v.empty())
      return 0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
  }
}

int main(int argc, char **argv)
{
  const std::string dir = arg(argc, argv, "--dir", ".");
  const std::string ch = arg(argc, argv, "--chunk", "0");
  const int nb = atoi(arg(argc, argv, "--nb", "36"));
  const int P = atoi(arg(argc, argv, "--P", "945"));
  const int N = atoi(arg(argc, argv, "--N", "3406"));
  const int L = atoi(arg(argc, argv, "--L", "28"));
  const int D = atoi(arg(argc, argv, "--digits", std::to_string(4 * L).c_str()));
  const int reps = atoi(arg(argc, argv, "--reps", "6"));
  const int preps = atoi(arg(argc, argv, "--profile-reps", "3"));
  const int n = (P + 31) / 32 * 32;
  Factor_Options fo;
  fo.equilibrate = atoi(arg(argc, argv, "--equilibrate", "1"));
  fo.panel_refine = atoi(arg(argc, argv, "--panel-refine", "1"));
  fo.refine_min_n = atoi(arg(argc, argv, "--refine-min", "0"));
  fo.correction_digits = atoi(arg(argc, argv, "--correction-digits", "0"));
  fo.symmetric_update_min = atoi(arg(argc, argv, "--gram-min", "0"));
  fo.skip_triangular_zeros = atoi(arg(argc, argv, "--skip-triangular-zeros", "1"));
  fo.lower_update = atoi(arg(argc, argv, "--lower-update", "1"));
  fo.correction_limit = atoi(arg(argc, argv, "--correction-limit", "0"));
  fo.structured_solve = atoi(arg(argc, argv, "--structured-solve", "1"));
  fo.reuse_residues = atoi(arg(argc, argv, "--reuse-residues", "1"));
  fo.refine_max_n = atoi(arg(argc, argv, "--refine-max", "0"));
  const std::string tag = arg(argc, argv, "--tag", "run");

  const int64_t sp = int64_t(nb) * P * P, bp = int64_t(nb) * P * N;
  auto Sm = slurp<int32_t>(dir + "/S" + ch + ".m.bin", size_t(L) * sp);
  auto Se = slurp<int64_t>(dir + "/S" + ch + ".e.bin", sp);
  auto Bm = slurp<int32_t>(dir + "/B" + ch + ".m.bin", size_t(L) * bp);
  auto Be = slurp<int64_t>(dir + "/B" + ch + ".e.bin", bp);

  // the structure of B: exact zero entries and zero columns, per block
  {
    int64_t zero = 0, zcols = 0;
    std::vector<int64_t> firstnz;
    for(int b = 0; b < nb; ++b)
      for(int j = 0; j < N; ++j)
        {
          int first = -1;
          for(int i = 0; i < P; ++i)
            if(Be[int64_t(b) * P * N + int64_t(i) * N + j] == MPF_EZ)
              ++zero;
            else if(first < 0)
              first = i;
          if(first < 0)
            ++zcols;
          else
            firstnz.push_back(first);
        }
    std::sort(firstnz.begin(), firstnz.end());
    printf("BSTRUCT %s: zero entries %.4f, zero columns %.4f of %ld; first "
           "nonzero row of a nonzero column, quartiles %ld %ld %ld of %d\n",
           tag.c_str(), double(zero) / bp, double(zcols) / (int64_t(nb) * N),
           long(int64_t(nb) * N),
           firstnz.empty() ? 0 : long(firstnz[firstnz.size() / 4]),
           firstnz.empty() ? 0 : long(firstnz[firstnz.size() / 2]),
           firstnz.empty() ? 0 : long(firstnz[3 * firstnz.size() / 4]), P);
  }

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("device %s; %s %s: %d blocks of %d (padded to %d), N=%d\n", prop.name,
         tag.c_str(), dir.c_str(), nb, P, n, N);

  // the stream-ordered pool keeps what is freed instead of returning it to
  // the device at every synchronization (as --retain-pool does)
  if(atoi(arg(argc, argv, "--retain-pool", "0")))
    {
      cudaMemPool_t pool;
      cudaDeviceGetDefaultMemPool(&pool, 0);
      uint64_t bytes = 64ULL << 30;
      cudaMemPoolSetAttribute(pool, cudaMemPoolAttrReleaseThreshold, &bytes);
    }
  Gram_Options base;
  base.prime_group = 32;
  base.fused_residues = true;
  base.crt_guard_bits = atoi(arg(argc, argv, "--crt-guard", "0"));
  base.share_workspace = atoi(arg(argc, argv, "--share-workspace", "1"));
  Block_Factor G(n, nb, L, D, base);
  // S: row-major P x P, uploaded as column-major (S symmetric); B: the dump's
  // row-major P x N is a column-major N x P, transposed. Padded with the
  // identity (S) and zero rows (B), as real_block_check and q_check.
  Mpf_Batch S0 = G.alloc(P, P), Bt = G.alloc(N, P), B0 = G.alloc(P, N);
  cudaMemcpy(S0.m, Sm.data(), Sm.size() * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(S0.e, Se.data(), Se.size() * 8, cudaMemcpyHostToDevice);
  cudaMemcpy(Bt.m, Bm.data(), Bm.size() * 4, cudaMemcpyHostToDevice);
  cudaMemcpy(Bt.e, Be.data(), Be.size() * 8, cudaMemcpyHostToDevice);
  mpf_transpose(Bt, B0, 0);
  G.release(Bt);
  Mpf_Batch S = G.alloc(n, n), B = G.alloc(n, N);
  {
    std::vector<int32_t> hm(size_t(L) * S.entries(), 0);
    std::vector<int64_t> he(S.entries(), MPF_EZ);
    for(int b = 0; b < nb; ++b)
      for(int i = 0; i < n; ++i)
        {
          const int64_t idx = int64_t(b) * n * n + int64_t(i) * n + i;
          hm[size_t(L - 1) * S.entries() + idx] = 1 << 29;
          he[idx] = -(32LL * L - 3);
        }
    cudaMemcpy(S.m, hm.data(), hm.size() * 4, cudaMemcpyHostToDevice);
    cudaMemcpy(S.e, he.data(), he.size() * 8, cudaMemcpyHostToDevice);
  }
  mpf_sub_block(S, 0, 0, S0, true, 0);
  mpf_zero(B, 0);
  mpf_sub_block(B, 0, 0, B0, true, 0);
  G.release(S0);
  G.release(B0);
  std::vector<int32_t>().swap(Sm);
  std::vector<int64_t>().swap(Se);
  std::vector<int32_t>().swap(Bm);
  std::vector<int64_t>().swap(Be);
  if(cudaDeviceSynchronize() != cudaSuccess)
    return printf("CUDA error after upload\n"), 2;
  Mpf_Batch Lf = G.alloc(n, n), X = G.alloc(n, n), Pn;

  // the diagnostic panel log first, once: it is not a timing
  if(atoi(arg(argc, argv, "--log-panels", "0")))
    {
      Block_Factor F(n, nb, L, D, base);
      Factor_Options o = fo;
      o.log_panels = true;
      o.measure_residual = atoi(arg(argc, argv, "--measure-residual", "0"));
      printf("PANELLOG begin %s\n", tag.c_str());
      F.factor(S, Lf, X, o);
      printf("PANELLOG end %s\n", tag.c_str());
      cudaDeviceSynchronize();
    }

  const auto run = [&](const char *label, Gram_Options opt, int count) {
    Block_Factor F(n, nb, L, D, opt);
    std::vector<double> tf, tp;
    double factor_high = 0, solve_high = 0;
    for(int rep = 0; rep < count; ++rep)
      {
        if(rep == 1)
          F.reset_timings();
        cudaDeviceSynchronize();
        pool_high();
        auto t0 = std::chrono::steady_clock::now();
        F.factor(S, Lf, X, fo);
        cudaDeviceSynchronize();
        auto t1 = std::chrono::steady_clock::now();
        factor_high = pool_high();
        if(Pn.m)
          G.release(Pn);
        F.solve(Lf, X, B, Pn, 0);
        cudaDeviceSynchronize();
        auto t2 = std::chrono::steady_clock::now();
        solve_high = pool_high();
        printf("%s %s rep %d: factor %.4f s, X B %.4f s; panels uncorrected %d, "
               "restarts %d\n",
               label, tag.c_str(), rep, seconds(t0, t1), seconds(t1, t2),
               F.panels_skipped(), F.correction_restarts());
        if(rep)
          {
            tf.push_back(seconds(t0, t1));
            tp.push_back(seconds(t1, t2));
          }
      }
    F.check();
    printf("%sMEDIAN %s factor %.4f X_B %.4f total %.4f\n", label, tag.c_str(),
           median(tf), median(tp), median(tf) + median(tp));
    for(const auto &kv : F.timings())
      printf("%s %s ms/rep %.3f : %s\n", label, tag.c_str(),
             kv.second / std::max(1, count - 1), kv.first.c_str());
    size_t free_b = 0, total_b = 0;
    cudaDeviceSynchronize();
    cudaMemGetInfo(&free_b, &total_b);
    printf("EVICTIONS %s %s %zu; products hold %.2f GB; device free %.2f of %.2f GB\n",
           label, tag.c_str(), F.cache_evictions(), F.product_bytes() / 1e9,
           free_b / 1e9, total_b / 1e9);
    // peak device use by stage: the pool's high-water mark (S, B, L, X, P and
    // temporaries) plus the products' cudaMalloc buffers
    printf("MEMORY %s %s device GB: factor %.2f, X B %.2f (pool high %.2f, %.2f; "
           "products %.2f); pool reserved high %.2f; host peak %.2f GB\n",
           label, tag.c_str(), (factor_high + F.product_bytes()) / 1e9,
           (solve_high + F.product_bytes()) / 1e9, factor_high / 1e9,
           solve_high / 1e9, F.product_bytes() / 1e9, pool_reserved_high() / 1e9,
           host_high() / 1e9);
    F.clear_cache();
    if(cudaDeviceSynchronize() != cudaSuccess)
      throw std::runtime_error("CUDA error");
  };

  // FNV-1a over the words of L, X and P (zero exponents included), to compare
  // builds bit for bit
  const auto digest = [&](const char *what, const Mpf_Batch &M) {
    std::vector<int32_t> hm(size_t(L) * M.entries());
    std::vector<int64_t> he(M.entries());
    cudaMemcpy(hm.data(), M.m, hm.size() * 4, cudaMemcpyDeviceToHost);
    cudaMemcpy(he.data(), M.e, he.size() * 8, cudaMemcpyDeviceToHost);
    unsigned long long h = 0xcbf29ce484222325ULL;
    for(int32_t v : hm)
      h = (h ^ (unsigned long long)(uint32_t)v) * 0x100000001b3ULL;
    for(int64_t v : he)
      h = (h ^ (unsigned long long)v) * 0x100000001b3ULL;
    printf("HASH %s %s %016llx\n", tag.c_str(), what, h);
  };

  // solve_flat against mm_flat, entry by entry, on the unpadded shape the solver uses
  const auto flat_check = [&]() {
    Block_Factor F(n, nb, L, D, base);
    F.factor(S, Lf, X, fo);
    Mpf_Batch Xc = G.alloc(P, P), Xt = G.alloc(P, P), Bc = G.alloc(P, N);
    mpf_sub_block(X, 0, 0, Xc, false, 0);
    mpf_transpose(Xc, Xt, 0);
    mpf_sub_block(B, 0, 0, Bc, false, 0);
    Mpf_Batch A = G.alloc(P, N), C = G.alloc(P, N);
    F.mm_flat(Xt, Bc, A.m, A.e);
    F.check();
    cudaDeviceSynchronize();
    auto t0 = std::chrono::steady_clock::now();
    F.solve_flat(Xt, Bc, C.m, C.e);
    F.check();
    cudaDeviceSynchronize();
    const double ts = seconds(t0, std::chrono::steady_clock::now());
    std::vector<int32_t> am(size_t(L) * A.entries()), cm(am.size());
    std::vector<int64_t> ae(A.entries()), ce(ae.size());
    cudaMemcpy(am.data(), A.m, am.size() * 4, cudaMemcpyDeviceToHost);
    cudaMemcpy(cm.data(), C.m, cm.size() * 4, cudaMemcpyDeviceToHost);
    cudaMemcpy(ae.data(), A.e, ae.size() * 8, cudaMemcpyDeviceToHost);
    cudaMemcpy(ce.data(), C.e, ce.size() * 8, cudaMemcpyDeviceToHost);
    size_t diff = 0, zeros = 0;
    for(size_t q = 0; q < ae.size(); ++q)
      {
        bool same = ae[q] == ce[q];
        for(int t = 0; t < L && same; ++t)
          same = am[size_t(t) * ae.size() + q] == cm[size_t(t) * ae.size() + q];
        diff += !same;
        zeros += ae[q] == MPF_EZ;
      }
    printf("FLATCHECK %s: solve_flat %.4f s; %zu of %zu entries differ from "
           "mm_flat; exact zeros %.4f\n",
           tag.c_str(), ts, diff, ae.size(), double(zeros) / ae.size());
    for(Mpf_Batch *M : {&Xc, &Xt, &Bc, &A, &C})
      G.release(*M);
  };

  auto stages = base;
  stages.timing = true;
  auto kernels = stages;
  kernels.profile_products = true;
  run("TIME", base, reps);
  if(atoi(arg(argc, argv, "--hash", "0")))
    {
      digest("L", Lf);
      digest("X", X);
      digest("P", Pn);
    }
  run("STAGE", stages, preps);
  run("KERNEL", kernels, preps);
  if(atoi(arg(argc, argv, "--flat-check", "0")))
    flat_check();
  run("TIME2", base, reps); // the same again: has the machine stayed the same?
  printf("PROFILE DONE %s\n", tag.c_str());
  return 0;
}
