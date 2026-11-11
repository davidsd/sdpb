// Distributed-Q collectives. Solver-facing header: no CUDA, NCCL or MPI types.
// NCCL is loaded only for G > 1; one process owns one device. SDPB supplies
// the small host broadcast used to exchange NCCL's unique ID.
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>

namespace sdpb_gpu
{
  struct Gpu_World
  {
    int rank = 0, size = 1, device = 0;
    std::function<void(void *, size_t, int)> host_broadcast;
  };

  class Gpu_Collective
  {
  public:
    explicit Gpu_Collective(const Gpu_World &world);
    ~Gpu_Collective();
    Gpu_Collective(const Gpu_Collective &) = delete;
    Gpu_Collective &operator=(const Gpu_Collective &) = delete;
    int rank() const;
    int size() const;
    // All device operations use stream 0, ordered with the backend's kernels.
    void broadcast(void *buffer, size_t bytes, int root);
    void max_u64(uint64_t *values, size_t count);
    void reduce_scatter_u32(const uint32_t *send, uint32_t *recv,
                            size_t recv_count);
    void allgather(const void *send, void *recv, size_t bytes);
    // Logical payload received, not a claim about physical link traffic.
    uint64_t received_bytes() const;
    uint64_t collective_calls() const;
  private:
    struct Impl;
    std::unique_ptr<Impl> p_;
  };
}
