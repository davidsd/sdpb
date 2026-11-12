#include "Gpu_Collective.hxx"
#include <cuda_runtime.h>
#include <dlfcn.h>
#include <unistd.h>
#include <chrono>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstring>
#include <array>

#ifdef SDPB_GPU_NCCL
#include <nccl.h>
#endif

namespace sdpb_gpu
{
  struct Gpu_Collective::Impl
  {
    Gpu_World w;
    uint64_t bytes = 0, calls = 0;
    bool timing = false;
    double t_broadcast = 0, t_max = 0, t_scatter = 0, t_gather = 0;
    using Clock = std::chrono::steady_clock;
    Clock::time_point begin()
    {
      if(!timing) return {};
      auto e = cudaStreamSynchronize(0);
      if(e != cudaSuccess) throw std::runtime_error(cudaGetErrorString(e));
      return Clock::now();
    }
    void end(double &total, Clock::time_point t)
    {
      if(!timing) return;
      auto e = cudaStreamSynchronize(0);
      if(e != cudaSuccess) throw std::runtime_error(cudaGetErrorString(e));
      total += std::chrono::duration<double>(Clock::now() - t).count();
    }
#ifdef SDPB_GPU_NCCL
    void *lib = nullptr;
    ncclComm_t comm = nullptr;
    decltype(&ncclGetUniqueId) get_id;
    decltype(&ncclCommInitRank) init;
    decltype(&ncclCommDestroy) destroy;
    decltype(&ncclBroadcast) bcast;
    decltype(&ncclAllReduce) allreduce;
    decltype(&ncclReduceScatter) scatter;
    decltype(&ncclAllGather) gather;
    decltype(&ncclGetErrorString) error;
    void check(ncclResult_t r)
    {
      if(r != ncclSuccess)
        throw std::runtime_error(std::string("NCCL: ") + error(r));
    }
    template<class T> void symbol(T &fn, const char *name)
    {
      fn = reinterpret_cast<T>(dlsym(lib, name));
      if(!fn) throw std::runtime_error(std::string("NCCL missing ") + name);
    }
#endif
    ~Impl()
    {
#ifdef SDPB_GPU_NCCL
      if(comm) destroy(comm);
      if(lib) dlclose(lib);
#endif
    }
  };

  Gpu_Collective::Gpu_Collective(const Gpu_World &w) : p_(new Impl)
  {
    if(w.size < 1 || w.rank < 0 || w.rank >= w.size)
      throw std::invalid_argument("Gpu_World rank/size");
    p_->w = w;
    const char *timing=std::getenv("SDPB_GPU_COMM_TIMING");
    p_->timing=timing && std::atoi(timing)!=0;
    int device=w.device;
    if(device<0)
      {
        int count=0;
        auto err=cudaGetDeviceCount(&count);
        if(err!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(err));
        // This rank's place among the ranks on its node, by host name.
        int local=w.local_rank>=0 ? w.local_rank : w.rank;
        if(w.size>1 && w.local_rank<0)
          {
            if(!w.host_broadcast)
              throw std::invalid_argument("Gpu_World needs a host broadcast");
            std::array<char,256> mine{};
            gethostname(mine.data(),mine.size()-1);
            local=0;
            for(int root=0;root<w.size;++root)
              {
                auto name=mine;
                w.host_broadcast(name.data(),name.size(),root);
                if(root<w.rank && name==mine) ++local;
              }
          }
        if(count!=1 && count<=local)
          throw std::runtime_error("fewer visible GPUs than MPI ranks on a node");
        device=count==1 ? 0 : local;
      }
    p_->w.device=device;
    auto e = cudaSetDevice(device);
    if(e != cudaSuccess) throw std::runtime_error(cudaGetErrorString(e));
    if(w.size == 1) return;
    if(!w.host_broadcast)
      throw std::invalid_argument("Gpu_World needs a host broadcast");
    // Slurm may expose all allocated GPUs, or one different GPU to each rank.
    // Refuse duplicate devices before NCCL can hang in communicator setup.
    cudaDeviceProp props;
    e=cudaGetDeviceProperties(&props,device);
    if(e!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(e));
    std::vector<cudaUUID_t> ids(w.size);
    for(int root=0;root<w.size;++root)
      {
        ids[root]=props.uuid;
        w.host_broadcast(&ids[root],sizeof(cudaUUID_t),root);
        for(int prev=0;prev<root;++prev)
          if(!std::memcmp(&ids[root],&ids[prev],sizeof(cudaUUID_t)))
            throw std::runtime_error("two MPI ranks selected the same GPU");
      }
#ifdef SDPB_GPU_NCCL
    const char *path = std::getenv("SDPB_GPU_NCCL_LIBRARY");
    p_->lib = dlopen(path ? path : "libnccl.so.2", RTLD_NOW | RTLD_LOCAL);
    if(!p_->lib) throw std::runtime_error(std::string("load NCCL: ") + dlerror());
    p_->symbol(p_->error, "ncclGetErrorString");
    p_->symbol(p_->get_id, "ncclGetUniqueId");
    p_->symbol(p_->init, "ncclCommInitRank");
    p_->symbol(p_->destroy, "ncclCommDestroy");
    p_->symbol(p_->bcast, "ncclBroadcast");
    p_->symbol(p_->allreduce, "ncclAllReduce");
    p_->symbol(p_->scatter, "ncclReduceScatter");
    p_->symbol(p_->gather, "ncclAllGather");
    ncclUniqueId id;
    if(w.rank == 0) p_->check(p_->get_id(&id));
    w.host_broadcast(&id, sizeof(id), 0);
    p_->check(p_->init(&p_->comm, w.size, id, w.rank));
#else
    throw std::runtime_error("multiple GPUs require a build with NCCL headers");
#endif
  }
  Gpu_Collective::~Gpu_Collective() = default;
  int Gpu_Collective::rank() const { return p_->w.rank; }
  int Gpu_Collective::size() const { return p_->w.size; }
  uint64_t Gpu_Collective::received_bytes() const { return p_->bytes; }
  uint64_t Gpu_Collective::collective_calls() const { return p_->calls; }
  bool Gpu_Collective::timing() const { return p_->timing; }
  double Gpu_Collective::seconds() const
  { return p_->t_broadcast + p_->t_max + p_->t_scatter + p_->t_gather; }
  double Gpu_Collective::broadcast_seconds() const { return p_->t_broadcast; }
  double Gpu_Collective::max_seconds() const { return p_->t_max; }
  double Gpu_Collective::reduce_scatter_seconds() const { return p_->t_scatter; }
  double Gpu_Collective::allgather_seconds() const { return p_->t_gather; }

  void Gpu_Collective::broadcast(void *buf, size_t bytes, int root)
  {
    if(root < 0 || root >= size()) throw std::invalid_argument("broadcast root");
    if(size() == 1 || !bytes) return;
#ifdef SDPB_GPU_NCCL
    const auto t = p_->begin();
    p_->check(p_->bcast(buf, buf, bytes, ncclUint8, root, p_->comm, 0));
    p_->end(p_->t_broadcast, t);
#endif
    ++p_->calls;
    if(rank() != root) p_->bytes += bytes;
  }
  void Gpu_Collective::max_u64(uint64_t *v, size_t n)
  {
    if(size() == 1 || !n) return;
#ifdef SDPB_GPU_NCCL
    const auto t = p_->begin();
    p_->check(p_->allreduce(v, v, n, ncclUint64, ncclMax, p_->comm, 0));
    p_->end(p_->t_max, t);
#endif
    ++p_->calls;
    p_->bytes += n * 8; // resulting payload, implementation-dependent wire volume
  }
  void Gpu_Collective::reduce_scatter_u32(const uint32_t *in, uint32_t *out,
                                         size_t n)
  {
    if(size() == 1)
      {
        auto e = cudaMemcpyAsync(out, in, n * 4, cudaMemcpyDeviceToDevice, 0);
        if(e != cudaSuccess) throw std::runtime_error(cudaGetErrorString(e));
        return;
      }
#ifdef SDPB_GPU_NCCL
    const auto t = p_->begin();
    p_->check(p_->scatter(in, out, n, ncclUint32, ncclSum, p_->comm, 0));
    p_->end(p_->t_scatter, t);
#endif
    ++p_->calls;
    p_->bytes += n * 4;
  }
  void Gpu_Collective::allgather(const void *in, void *out, size_t n)
  {
    if(size() == 1)
      {
        auto e = cudaMemcpyAsync(out, in, n, cudaMemcpyDeviceToDevice, 0);
        if(e != cudaSuccess) throw std::runtime_error(cudaGetErrorString(e));
        return;
      }
#ifdef SDPB_GPU_NCCL
    const auto t = p_->begin();
    p_->check(p_->gather(in, out, n, ncclUint8, p_->comm, 0));
    p_->end(p_->t_gather, t);
#endif
    ++p_->calls;
    p_->bytes += n * (size() - 1);
  }
}
