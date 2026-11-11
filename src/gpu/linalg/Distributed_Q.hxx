// Distributed exact Q and corrected tiled Cholesky.
// Uses SDPB's Gram/Block_Factor primitives, based on Antoine Tilloy's SDPviBe;
// see docs/GPU_backend_credits.md. No whole Q or whole inverse is replicated.
#pragma once
#include "Block_Factor.hxx"
#include "../solver/Gpu_Collective.hxx"
#include <memory>
#include <vector>
#include <string>

namespace sdpb_gpu
{
  struct Distributed_Q_Options
  {
    int tile = 512;
    int row_panel = 8192; // bounds residue planes independently of chunk size
    int refine = 1;
    bool wide_residual = false; // experiment: one rounded integer reduction
  };
  class Distributed_Q
  {
  public:
    Distributed_Q(int N, int L, int digits, int64_t total_rows,
                   Gpu_Collective &world, Gram_Options gram_options = {},
                   Distributed_Q_Options options = {});
    ~Distributed_Q();
    // All ranks call with their local P chunks in flat column-major layout.
    void assemble(const std::vector<Mpf_Matrix> &local_P);
    // Test seam: a full host Q is read in tiles, by each tile's owner only.
    void load_tile(int i, int j, const int32_t *m, const int64_t *e);
    void download_tile(int i, int j, int32_t *m, int64_t *e, bool factor = false);
    void factor();
    void solve(const Mpf_Batch &rhs, Mpf_Batch &out, int refine = -1);
    void multiply(const Mpf_Batch &v, Mpf_Batch &out);
    // rhs-Q*v: full products, guarded integer rank sum, one output rounding.
    void residual(const Mpf_Batch &rhs, const Mpf_Batch &v, Mpf_Batch &out);
    void forward(const Mpf_Batch &rhs, Mpf_Batch &out);
    void factor_diagonal(Mpf_Batch &out);
    // Per-tile byte digests, independent of rank count. No all-gather of Q.
    std::string digests(bool factor = false);
    void save_factor(const std::string &prefix); // test export: owner tiles
    void corrupt_tile(int i, int j, int bit); // negative control of factor
    // Graded control: clears the low limbs of every entry of a factor tile
    void corrupt_tile_limbs(int i, int j, int limbs);
    int tile_size() const;
    int tiles() const;
    int owner(int i, int j) const;
    uint64_t stored_bytes() const;
    uint64_t peak_device_bytes() const;
    // Inner phase timing resets CUDA's high-water counter. Forward this peak
    // to the outer assembly/factor/solve stage before it closes its own epoch.
    uint64_t last_operation_peak_device_bytes() const;
    std::string profile_json() const;
  private:
    struct Impl;
    std::unique_ptr<Impl> p_;
  };

  // Sum device vectors, same order on all ranks; O(G*N*E) transient memory.
  // This is floating arithmetic, unlike the exact modular Q reduction.
  void distributed_sum(const Mpf_Batch &local, Mpf_Batch &sum,
                        Gpu_Collective &world);
  // Optional validation of the unchanged G=1 path against distributed tiles.
  // Uses one tile of staging; the full Q never moves to host memory.
  std::string full_Q_tile_digests(const Mpf_Batch &Q,int tile);
}
