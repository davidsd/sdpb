// NCCL and host MPI rates between this job's ranks, at
// the sizes distributed Q exchanges. One process per GPU; ranks may share a
// node. Each RATE line: operation, bytes per rank, repetitions, best and median
// seconds, and bytes/best in GB/s (decimal). Point-to-point, reduce-scatter,
// broadcast and all-gather check their data after the last repetition, so a rate
// from a transport that corrupts data can't pass; all-reduce is only timed.
#include <cuda_runtime.h>
#include <dlfcn.h>
#include <mpi.h>
#include <nccl.h>
#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#define CK(x) do { auto e_=(x); if(e_!=cudaSuccess) throw std::runtime_error(std::string(#x)+": "+cudaGetErrorString(e_)); } while(0)

namespace
{
  decltype(&ncclCommInitRank) f_init; decltype(&ncclGetUniqueId) f_id;
  decltype(&ncclSend) f_send; decltype(&ncclRecv) f_recv;
  decltype(&ncclBroadcast) f_bcast; decltype(&ncclAllReduce) f_allreduce;
  decltype(&ncclReduceScatter) f_rs; decltype(&ncclAllGather) f_gather;
  decltype(&ncclGroupStart) f_gs; decltype(&ncclGroupEnd) f_ge;
  decltype(&ncclGetErrorString) f_err; decltype(&ncclGetVersion) f_ver;
  decltype(&ncclCommDestroy) f_destroy;
  void N(ncclResult_t r) { if(r!=ncclSuccess) throw std::runtime_error(std::string("NCCL: ")+f_err(r)); }
  template<class T> void sym(void *l,T &f,const char *n)
  { f=reinterpret_cast<T>(dlsym(l,n)); if(!f) throw std::runtime_error(std::string("missing ")+n); }

  int rank_=0,size_=1;
  double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

  // Times op over reps (after one warm-up); the slowest rank's time counts.
  template<class F> void rate(const char *name,size_t bytes,int reps,F op,const char *extra="")
  {
    op(); CK(cudaStreamSynchronize(0)); MPI_Barrier(MPI_COMM_WORLD);
    std::vector<double> t;
    for(int r=0;r<reps;++r)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        const double t0=now(); op(); CK(cudaStreamSynchronize(0));
        double dt=now()-t0,mx=0;
        MPI_Allreduce(&dt,&mx,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        t.push_back(mx);
      }
    std::sort(t.begin(),t.end());
    if(!rank_) { printf("RATE %-22s bytes=%-11zu reps=%d best=%.6f med=%.6f GBps=%.3f %s\n",name,bytes,reps,t[0],t[t.size()/2],bytes/t[0]/1e9,extra); fflush(stdout); }
  }
  // Fills a device buffer with value+index, as 32-bit words.
  void fill(uint32_t *d,size_t words,uint32_t value,std::vector<uint32_t> &h)
  {
    h.resize(words); for(size_t i=0;i<words;++i) h[i]=value+uint32_t(i);
    CK(cudaMemcpy(d,h.data(),words*4,cudaMemcpyHostToDevice));
  }
  void expect(const uint32_t *d,size_t words,const std::vector<uint32_t> &want,const char *what)
  {
    std::vector<uint32_t> h(words); CK(cudaMemcpy(h.data(),d,words*4,cudaMemcpyDeviceToHost));
    if(std::memcmp(h.data(),want.data(),words*4))
      throw std::runtime_error(std::string("data differs after ")+what);
    printf("CHECK rank=%d %s bytes=%zu ok\n",rank_,what,words*4); fflush(stdout);
  }
}

int main(int argc,char **argv)
{
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_); MPI_Comm_size(MPI_COMM_WORLD,&size_);
  try
    {
      const double max_gb=argc>1 ? atof(argv[1]) : 2.1;
      char host[256]={}; gethostname(host,255);
      std::vector<char> hosts(256*size_);
      MPI_Allgather(host,256,MPI_CHAR,hosts.data(),256,MPI_CHAR,MPI_COMM_WORLD);
      std::map<std::string,int> on_node; int local=0;
      for(int r=0;r<size_;++r)
        {
          const std::string h(&hosts[256*r]);
          if(r==rank_) local=on_node[h];
          on_node[h]++;
        }
      int devices=0; CK(cudaGetDeviceCount(&devices));
      const int device=devices==1?0:local;
      CK(cudaSetDevice(device));
      cudaDeviceProp prop; CK(cudaGetDeviceProperties(&prop,device));
      printf("RANK %d host=%s local=%d device=%d of %d name=%s\n",rank_,host,local,device,devices,prop.name); fflush(stdout);
      const char *path=getenv("SDPB_GPU_NCCL_LIBRARY");
      void *lib=dlopen(path?path:"libnccl.so.2",RTLD_NOW|RTLD_LOCAL);
      if(!lib) throw std::runtime_error(dlerror());
      sym(lib,f_init,"ncclCommInitRank"); sym(lib,f_id,"ncclGetUniqueId"); sym(lib,f_send,"ncclSend"); sym(lib,f_recv,"ncclRecv");
      sym(lib,f_bcast,"ncclBroadcast"); sym(lib,f_allreduce,"ncclAllReduce"); sym(lib,f_rs,"ncclReduceScatter"); sym(lib,f_gather,"ncclAllGather");
      sym(lib,f_gs,"ncclGroupStart"); sym(lib,f_ge,"ncclGroupEnd"); sym(lib,f_err,"ncclGetErrorString"); sym(lib,f_ver,"ncclGetVersion"); sym(lib,f_destroy,"ncclCommDestroy");
      int ver=0; f_ver(&ver); if(!rank_) printf("NCCL version %d\n",ver);
      ncclUniqueId id; if(!rank_) N(f_id(&id));
      MPI_Bcast(&id,sizeof(id),MPI_BYTE,0,MPI_COMM_WORLD);
      ncclComm_t comm;
      double t0=now(); N(f_init(&comm,size_,id,rank_)); if(!rank_) { printf("INIT %.3f s\n",now()-t0); fflush(stdout); }
      const size_t top=size_t(max_gb*1e9/4)*4;
      void *dbuf=nullptr,*dsend=nullptr;
      CK(cudaMalloc(&dbuf,top+4096)); CK(cudaMalloc(&dsend,(top+4096)*size_));
      auto *a=static_cast<uint32_t*>(dbuf); auto *s=static_cast<uint32_t*>(dsend);
      std::vector<uint32_t> h,want;
      std::vector<size_t> sizes={4096,65536,size_t(1)<<20,size_t(8)<<20,size_t(32)<<20,size_t(128)<<20,size_t(512)<<20,size_t(2048)<<20};
      sizes.erase(std::remove_if(sizes.begin(),sizes.end(),[&](size_t b){return b>top;}),sizes.end());
      auto reps=[](size_t b){ return b<=(size_t(8)<<20) ? 20 : b<=(size_t(128)<<20) ? 6 : 3; };
      // Point to point, rank 0 to each other rank, then both ways together.
      for(int peer=1;peer<size_;++peer)
        {
          const std::string tag=!strcmp(&hosts[0],&hosts[256*peer])?"intra":"inter";
          const std::string nm="p2p-"+tag+"-0to"+std::to_string(peer),nb="bidir-"+tag+"-0to"+std::to_string(peer);
          const bool in=rank_==0||rank_==peer;
          for(size_t b:sizes)
            {
              MPI_Barrier(MPI_COMM_WORLD);
              if(in) fill(s,b/4,7+peer,h);
              rate(nm.c_str(),b,reps(b),[&]{
                if(rank_==0) N(f_send(s,b,ncclUint8,peer,comm,0));
                else if(rank_==peer) N(f_recv(a,b,ncclUint8,0,comm,0));
              });
              if(rank_==peer) expect(a,b/4,h,nm.c_str());
              if(in) fill(s,b/4,rank_==0?100:200,h);
              rate(nb.c_str(),b,reps(b),[&]{
                if(!in) return;
                N(f_gs());
                N(f_send(s,b,ncclUint8,rank_==0?peer:0,comm,0));
                N(f_recv(a,b,ncclUint8,rank_==0?peer:0,comm,0));
                N(f_ge());
              });
              if(in) { std::vector<uint32_t> w(b/4); for(size_t i=0;i<b/4;++i) w[i]=(rank_==0?200:100)+uint32_t(i); expect(a,b/4,w,nb.c_str()); }
            }
        }
      // Collectives, with the sizes Q uses for each.
      for(size_t b:sizes)
        {
          const size_t w=b/4;
          // Each rank contributes value (rank+1)+index for each destination's segment.
          h.resize(w*size_);
          for(size_t i=0;i<h.size();++i) h[i]=uint32_t(rank_+1)+uint32_t(i%w);
          CK(cudaMemcpy(s,h.data(),h.size()*4,cudaMemcpyHostToDevice));
          rate("reduce_scatter_u32",b,reps(b),[&]{ N(f_rs(s,a,w,ncclUint32,ncclSum,comm,0)); });
          want.resize(w); for(size_t i=0;i<w;++i) want[i]=uint32_t(size_)*uint32_t(i)+uint32_t(size_*(size_+1)/2);
          expect(a,w,want,"reduce_scatter_u32");
          if(!rank_) fill(a,w,5,h);
          rate("broadcast",b,reps(b),[&]{ N(f_bcast(a,a,b,ncclUint8,0,comm,0)); });
          h.resize(w); for(size_t i=0;i<w;++i) h[i]=5+uint32_t(i);
          expect(a,w,h,"broadcast");
          if(size_*b<=top+4096)
            {
              const size_t each=b/size_ & ~size_t(3);
              fill(s,each/4,1000*(rank_+1),h);
              rate("allgather",each,reps(b),[&]{ N(f_gather(s,a,each,ncclUint8,comm,0)); });
              want.resize(each/4*size_);
              for(size_t i=0;i<want.size();++i) want[i]=1000*uint32_t(i/(each/4)+1)+uint32_t(i%(each/4));
              expect(a,want.size(),want,"allgather");
            }
          rate("allreduce_u64_max",b,reps(b),[&]{ N(f_allreduce(a,a,b/8,ncclUint64,ncclMax,comm,0)); });
        }
      // Host MPI, which carries the driver's scalars and host broadcasts.
      {
        std::vector<char> hb(size_t(256)<<20,1),hc(hb.size());
        for(size_t b:{size_t(8),size_t(4096),size_t(1)<<20,size_t(32)<<20,size_t(256)<<20})
          {
            rate("mpi_bcast",b,b<(size_t(1)<<20)?200:5,[&]{ MPI_Bcast(hb.data(),int(b),MPI_BYTE,0,MPI_COMM_WORLD); });
            if(size_>1) rate("mpi_sendrecv_0to1",b,b<(size_t(1)<<20)?200:5,[&]{
              if(rank_==0) MPI_Send(hb.data(),int(b),MPI_BYTE,1,0,MPI_COMM_WORLD);
              else if(rank_==1) MPI_Recv(hc.data(),int(b),MPI_BYTE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            });
          }
        double x=1,y=0;
        rate("mpi_allreduce_double",8,200,[&]{ MPI_Allreduce(&x,&y,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); });
        // Host and device, for the case where NCCL stages through host memory.
        void *pin=nullptr; size_t pb=size_t(1)<<30; CK(cudaMallocHost(&pin,pb));
        rate("host_to_device_pinned",pb,3,[&]{ CK(cudaMemcpy(a,pin,pb,cudaMemcpyHostToDevice)); });
        rate("device_to_host_pinned",pb,3,[&]{ CK(cudaMemcpy(pin,a,pb,cudaMemcpyDeviceToHost)); });
        CK(cudaFreeHost(pin));
      }
      CK(cudaFree(dbuf)); CK(cudaFree(dsend));
      N(f_destroy(comm));
      if(!rank_) { printf("DONE\n"); fflush(stdout); }
      MPI_Finalize(); return 0;
    }
  catch(const std::exception &e)
    { fprintf(stderr,"rank %d: %s\n",rank_,e.what()); fflush(stderr); MPI_Abort(MPI_COMM_WORLD,1); return 1; }
}
