// Tests for distributed Q: exact assembly against the existing full Gram,
// real-Q factor exports and forward solves for the independent MPFR gate.
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include "../linalg/Distributed_Q.hxx"
#include "../solver/Mpf_Ops.hxx"
#include "../rns/rns_kernels.cuh"
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>
#include <chrono>

using namespace sdpb_gpu;
#define CHECK(x) RNS_CUDA_CHECK(x)
constexpr int64_t ez=-(int64_t(1)<<40);
const char *arg(int argc,char **argv,const char *key,const char *def)
{ for(int i=1;i+1<argc;++i) if(!strcmp(argv[i],key)) return argv[i+1]; return def; }
struct Array
{
  Mpf_Batch b;
  Array(int L,int r,int c) : b{nullptr,nullptr,L,1,r,c}
  { CHECK(cudaMallocAsync(&b.m,size_t(L)*r*c*4,0)); CHECK(cudaMallocAsync(&b.e,size_t(r)*c*8,0)); }
  ~Array() { cudaFreeAsync(b.m,0); cudaFreeAsync(b.e,0); }
};
template<class T> struct Mapped
{
  const T *p=nullptr; size_t bytes=0;
  explicit Mapped(const std::string &file)
  {
    int fd=open(file.c_str(),O_RDONLY); struct stat st;
    if(fd<0 || fstat(fd,&st)) throw std::runtime_error("open "+file);
    bytes=st.st_size; p=static_cast<const T *>(mmap(nullptr,bytes,PROT_READ,MAP_PRIVATE,fd,0)); close(fd);
    if(p==MAP_FAILED) throw std::runtime_error("mmap "+file);
  }
  ~Mapped() { if(p && p!=MAP_FAILED) munmap(const_cast<T *>(p),bytes); }
};
void entry_one(std::vector<int32_t> &m,std::vector<int64_t> &e,int L,int64_t plane,int64_t q)
{ m[(L-1)*plane+q]=0x20000000; e[q]=-(32LL*L-3); }

void assembly(Gpu_Collective &world,int N,int tile,bool drop,int L)
{
  constexpr int nchunks=11;
  Gram_Options go; go.prime_group=32; go.fused_residues=true;
  go.max_rect=256; go.crt_batch=16384;
  int64_t total=0; for(int c=0;c<nchunks;++c) total+=37+3*c;
  Distributed_Q q(N,L,4*L,total,world,go,{tile,64,1});
  std::vector<std::unique_ptr<Array>> arrays;
  std::vector<Mpf_Matrix> local;
  std::unique_ptr<Mpf_Gram> full;
  if(world.rank()==0) full=std::make_unique<Mpf_Gram>(N,L,4*L,total,go);
  auto fill=[&](int c,Array &a) {
    const int64_t n=a.b.entries(); std::vector<int32_t> m(size_t(L)*n,0); std::vector<int64_t> e(n,ez);
    for(int j=0;j<N;++j) for(int i=0;i<a.b.rows;++i)
      {
        const int64_t k=int64_t(j)*a.b.rows+i;
        // Normalized, signed, nontrivial low bits and column scales.
        uint32_t x=uint32_t((c+13)*12345u+(j+3)*54321u+(i+7)*4321u);
        for(int t=0;t<L;++t) { x=1664525u*x+1013904223u; m[int64_t(t)*n+k]=int32_t(x); }
        m[(L-1)*n+k]=(uint32_t(m[(L-1)*n+k])&0x1fffffff)|0x20000000;
        e[k]=-(32LL*L-3)+(j%17-8)*20+(c%3)*3;
        if((c+i+j)%7==0) { for(int t=0;t<L;++t) m[int64_t(t)*n+k]=0; e[k]=ez; }
        else if((i+j)%3==0)
          { uint64_t carry=1; for(int t=0;t<L;++t) { uint64_t z=uint64_t(~uint32_t(m[int64_t(t)*n+k]))+carry; m[int64_t(t)*n+k]=uint32_t(z); carry=z>>32; } }
      }
    mpf_upload(m.data(),e.data(),a.b,0);
  };
  // Single-GPU reference receives all chunks, including the negative-control omission.
  std::vector<std::unique_ptr<Array>> reference;
  for(int c=0;c<nchunks;++c)
    {
      if(c%world.size()==world.rank() && !(drop && c==5))
        {
          arrays.emplace_back(new Array(L,37+3*c,N)); fill(c,*arrays.back());
          local.push_back(arrays.back()->b.view());
        }
      if(world.rank()==0)
        { reference.emplace_back(new Array(L,37+3*c,N)); fill(c,*reference.back()); full->scan(reference.back()->b.view()); }
    }
  q.assemble(local);
  const int64_t packed=int64_t(N)*(N+1)/2;
  std::vector<int32_t> rm(size_t(L)*packed);
  std::vector<int64_t> re(packed);
  if(world.rank()==0)
    {
      for(auto &a:reference) full->add(a->b.view());
      Array out(L,int(packed),1); full->finish(out.b.m,out.b.e);
      mpf_download(out.b,rm.data(),re.data(),0);
    }
  MPI_Bcast(rm.data(),int(rm.size()),MPI_INT32_T,0,MPI_COMM_WORLD);
  MPI_Bcast(re.data(),int(re.size()),MPI_INT64_T,0,MPI_COMM_WORLD);
  long long bad=0, checked=0;
  std::vector<int32_t> tm(size_t(L)*tile*tile); std::vector<int64_t> te(size_t(tile)*tile);
  for(int j=0;j<q.tiles();++j) for(int i=j;i<q.tiles();++i) if(q.owner(i,j)==world.rank())
    {
      q.download_tile(i,j,tm.data(),te.data());
      for(int col=0;col<tile;++col) for(int row=0;row<tile;++row)
        {
          int ii=i*tile+row,jj=j*tile+col; if(ii>=N || jj>=N) continue;
          const int64_t a=Gram_Accumulator::packed_index(N,std::max(ii,jj),std::min(ii,jj));
          const int64_t v=int64_t(col)*tile+row; bool wrong=te[v]!=re[a];
          for(int t=0;t<L;++t) wrong|=tm[int64_t(t)*tile*tile+v]!=rm[int64_t(t)*packed+a];
          bad+=wrong; ++checked;
        }
    }
  long long allbad=0, allchecked=0;
  MPI_Allreduce(&bad,&allbad,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&checked,&allchecked,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  if(world.rank()==0) std::cout << "ASSEMBLY N="<<N<<" G="<<world.size()<<" entries="<<allchecked<<" different="<<allbad<<"\n";
  std::cout << "PROFILE " << q.profile_json() << "\n";
  if(allbad) throw std::runtime_error("assembly differs from full single-GPU Gram");
}

template<class Solver>
void export_solutions(Solver &q,Gpu_Collective &world,int N,int columns,
                      const std::string &dir,const std::string &out,int L)
{
  Array rhs(L,N,1), result(L,N,1), solved(L,N,1);
  std::vector<int32_t> rm(size_t(L)*N), cm(size_t(L)*N);
  std::vector<int64_t> re(N), ce(N);
  Mapped<int32_t> bm(dir+"/B0.m.bin"), pm(dir+"/P0.m.bin");
  Mapped<int64_t> be(dir+"/B0.e.bin"), pe(dir+"/P0.e.bin");
  FILE *f=nullptr;
  if(world.rank()==0)
    {
      f=fopen((out+".cols.bin").c_str(),"wb"); if(!f) throw std::runtime_error("cols output");
      fwrite("QCOL",1,4,f); int32_t h[5]={1,N,32,L,columns}; fwrite(h,4,5,f);
      for(int j=0;j<columns;++j) { int32_t x=j; fwrite(&x,4,1,f); }
      for(int j=0;j<columns;++j) { int32_t x=0; fwrite(&x,4,1,f); }
    }
  for(int j=0;j<columns;++j)
    {
      for(int i=0;i<N;++i)
        { re[i]=be.p[size_t(i)*32+j]; for(int t=0;t<L;++t) rm[size_t(t)*N+i]=bm.p[(size_t(t)*N+i)*32+j]; }
      mpf_upload(rm.data(),re.data(),rhs.b,0);
      q.forward(rhs.b,result.b);
      // The full solve is exercised too; independent residual checks follow.
      q.solve(rhs.b,solved.b);
      mpf_download(result.b,cm.data(),ce.data(),0);
      if(f)
        {
          for(int which=0;which<3;++which) for(int i=0;i<N;++i)
            {
              for(int t=0;t<L;++t)
                { int32_t v=which==0?rm[size_t(t)*N+i]:(which==1?cm[size_t(t)*N+i]:pm.p[(size_t(t)*N+i)*32+j]); fwrite(&v,4,1,f); }
              int64_t e=which==0?re[i]:(which==1?ce[i]:pe.p[size_t(i)*32+j]); fwrite(&e,8,1,f);
            }
          std::vector<int32_t> sm(size_t(L)*N); std::vector<int64_t> se(N);
          mpf_download(solved.b,sm.data(),se.data(),0);
          std::ofstream fm(out+".solve"+std::to_string(j)+".m.bin",std::ios::binary), fe(out+".solve"+std::to_string(j)+".e.bin",std::ios::binary);
          fm.write(reinterpret_cast<char *>(sm.data()),sm.size()*4); fe.write(reinterpret_cast<char *>(se.data()),se.size()*8);
        }
    }
  if(f && fclose(f)) throw std::runtime_error("cols write");
}

void real_factor(Gpu_Collective &world,int argc,char **argv)
{
  const int N=atoi(arg(argc,argv,"--N","3406")), L=atoi(arg(argc,argv,"--L","28"));
  const int tile=atoi(arg(argc,argv,"--tile","512"));
  const int columns=atoi(arg(argc,argv,"--columns","32"));
  const bool corrupt=atoi(arg(argc,argv,"--corrupt","0"));
  std::string dir=arg(argc,argv,"--dir",""), out=arg(argc,argv,"--out","dq");
  if(columns<1 || columns>32) throw std::invalid_argument("columns must be 1..32 for frozen B0/P0");
  Gram_Options go; go.prime_group=32; go.fused_residues=true; go.defer_crt_check=true;
  go.crt_guard_bits=atoi(arg(argc,argv,"--crt-guard","0"));
  if(world.rank()==0) std::cout << "FROZEN_FACTOR crt_guard=" << go.crt_guard_bits << "\n";
  const bool wide=atoi(arg(argc,argv,"--wide-residual","0"));
  const int refinements=atoi(arg(argc,argv,"--refine","1"));
  if(world.rank()==0) std::cout << "RESIDUAL wide="<<wide<<" refine="<<refinements<<"\n";
  Distributed_Q q(N,L,4*L,1,world,go,{tile,8192,refinements,wide});
  {
    Mapped<int32_t> m(dir+"/S0.m.bin"); Mapped<int64_t> e(dir+"/S0.e.bin");
    if(m.bytes!=size_t(L)*N*N*4 || e.bytes!=size_t(N)*N*8) throw std::runtime_error("Q file shape");
    std::vector<int32_t> tm(size_t(L)*tile*tile); std::vector<int64_t> te(size_t(tile)*tile);
    for(int j=0;j<q.tiles();++j) for(int i=j;i<q.tiles();++i) if(q.owner(i,j)==world.rank())
      {
        std::fill(tm.begin(),tm.end(),0); std::fill(te.begin(),te.end(),ez);
        for(int col=0;col<tile;++col) for(int row=0;row<tile;++row)
          {
            int ii=i*tile+row,jj=j*tile+col; if(ii>=N || jj>=N) continue;
            size_t src=size_t(ii)*N+jj, dest=size_t(col)*tile+row;
            te[dest]=e.p[src]; for(int t=0;t<L;++t) tm[size_t(t)*tile*tile+dest]=m.p[size_t(t)*N*N+src];
          }
        q.load_tile(i,j,tm.data(),te.data());
      }
  }
  std::ofstream dg(out+".q.rank"+std::to_string(world.rank())+".digest"); dg<<q.digests();
  if(std::string(arg(argc,argv,"--mode","factor"))=="residual")
    {
      const std::string vp=arg(argc,argv,"--vector","");
      Mapped<int32_t> m(vp+".m.bin");Mapped<int64_t> e(vp+".e.bin");
      if(m.bytes!=size_t(L)*N*4 || e.bytes!=size_t(N)*8)throw std::runtime_error("residual vector shape");
      Array v(L,N,1),rhs(L,N,1),r(L,N,1);
      mpf_upload(m.p,e.p,v.b,0);
      std::vector<int32_t> rm(size_t(L)*N,0);std::vector<int64_t> re(N,ez);
      const int col=atoi(arg(argc,argv,"--rhs-index","0"));
      if(col<0 || col>=N)throw std::invalid_argument("rhs index");
      entry_one(rm,re,L,N,col);mpf_upload(rm.data(),re.data(),rhs.b,0);
      q.residual(rhs.b,v.b,r.b);
      mpf_download(r.b,rm.data(),re.data(),0);
      const std::string prefix=out+".residual.rank"+std::to_string(world.rank());
      std::ofstream fm(prefix+".m.bin",std::ios::binary),fe(prefix+".e.bin",std::ios::binary);
      fm.write(reinterpret_cast<char *>(rm.data()),rm.size()*4);
      fe.write(reinterpret_cast<char *>(re.data()),re.size()*8);
      if(!fm || !fe)throw std::runtime_error("residual export failed");
      return;
    }
  q.factor();
  if(corrupt) q.corrupt_tile(0,0,0);
  const int corrupt_limbs=atoi(arg(argc,argv,"--corrupt-limbs","0"));
  if(corrupt_limbs) q.corrupt_tile_limbs(0,0,corrupt_limbs);
  std::ofstream fdg(out+".L.rank"+std::to_string(world.rank())+".digest"); fdg<<q.digests(true);
  q.save_factor(out+".L");
  export_solutions(q,world,N,columns,dir,out,L);
  std::cout << "PROFILE " << q.profile_json() << "\n";
}

// Full-matrix factor and inverse-product solve, matching Gpu_Solver's G=1
// arithmetic, for a separate baseline of the full-Q solve error.
struct Monolithic_Q
{
  int N, L;
  Array Q, X;
  Monolithic_Q(int n,const std::string &dir,int guard,int L_) : N(n),L(L_),Q(L,n,n),X(L,n,n)
  {
    {
      Mapped<int32_t> m(dir+"/S0.m.bin"); Mapped<int64_t> e(dir+"/S0.e.bin");
      if(m.bytes!=size_t(L)*N*N*4 || e.bytes!=size_t(N)*N*8)
        throw std::runtime_error("monolithic Q shape");
      std::vector<int32_t> hm(size_t(L)*N*N);
      std::vector<int64_t> he(size_t(N)*N);
      for(int j=0;j<N;++j) for(int i=0;i<N;++i)
        {
          const size_t src=size_t(i)*N+j,dst=size_t(j)*N+i;
          he[dst]=e.p[src];
          for(int k=0;k<L;++k) hm[size_t(k)*N*N+dst]=m.p[size_t(k)*N*N+src];
        }
      mpf_upload(hm.data(),he.data(),Q.b,0);
    }
    const int np=(N+31)/32*32;
    Array padded(L,np,np),F(L,np,np),Xi(L,np,np);
    mpf_pad_identity(Q.b,padded.b,0);
    Gram_Options go; go.prime_group=32; go.fused_residues=true; go.defer_crt_check=true;
    go.crt_guard_bits=guard;
    std::cout << "MONOLITHIC_FACTOR crt_guard=" << guard << "\n";
    Block_Factor engine(32,1,L,4*L,go,0);
    Factor_Options fo{true,1,0,-1,1024,true};
    auto t=std::chrono::steady_clock::now();
    engine.factor(padded.b,F.b,Xi.b,fo);
    engine.check(); CHECK(cudaStreamSynchronize(0));
    std::cout << "MONOLITHIC_FACTOR seconds="
              << std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<"\n";
    mpf_sub_block(Xi.b,0,0,X.b,false,0);
  }
  void mv(const Mpf_Batch &A,bool trans,const Mpf_Batch &v,Mpf_Batch &out)
  {
    Gemv_Spec g; g.no=N; g.nti=N; g.ao=trans?N:1; g.ati=trans?1:N; g.vti=1;
    mpf_gemv(A,v,out,g,0);
  }
  void forward(const Mpf_Batch &rhs,Mpf_Batch &out) { mv(X.b,false,rhs,out); }
  void solve(const Mpf_Batch &rhs,Mpf_Batch &out)
  {
    Array temp(L,N,1),prod(L,N,1),residual(L,N,1),correction(L,N,1);
    mv(X.b,false,rhs,temp.b); mv(X.b,true,temp.b,out);
    mv(Q.b,false,out,prod.b); mpf_batch_add(rhs,prod.b,-1,residual.b,0);
    mv(X.b,false,residual.b,temp.b); mv(X.b,true,temp.b,correction.b);
    mpf_batch_add(out,correction.b,1,out,0);
  }
};
void monolithic(Gpu_Collective &world,int argc,char **argv)
{
  if(world.size()!=1) throw std::invalid_argument("monolithic baseline needs one rank");
  const int N=atoi(arg(argc,argv,"--N","3406"));
  const int columns=atoi(arg(argc,argv,"--columns","32"));
  if(N<=0 || columns<1 || columns>32) throw std::invalid_argument("monolithic dimensions");
  const std::string dir=arg(argc,argv,"--dir",""),out=arg(argc,argv,"--out","mono");
  const int L=atoi(arg(argc,argv,"--L","28"));
  Monolithic_Q q(N,dir,atoi(arg(argc,argv,"--crt-guard","0")),L);
  export_solutions(q,world,N,columns,dir,out,L);
}

int main(int argc,char **argv)
{
  MPI_Init(&argc,&argv);
  int rank=0,size=1; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
  try
    {
      int devices=0; CHECK(cudaGetDeviceCount(&devices));
      Gpu_World w; w.rank=rank; w.size=size; w.device=devices==1?0:rank%devices;
      w.host_broadcast=[](void *p,size_t n,int root){ MPI_Bcast(p,int(n),MPI_BYTE,root,MPI_COMM_WORLD); };
      {
        Gpu_Collective world(w);
        cudaMemPool_t pool; CHECK(cudaDeviceGetDefaultMemPool(&pool,w.device));
        uint64_t hold=UINT64_MAX; CHECK(cudaMemPoolSetAttribute(pool,cudaMemPoolAttrReleaseThreshold,&hold));
        const std::string mode=arg(argc,argv,"--mode","assembly");
        auto t=std::chrono::steady_clock::now();
        if(mode=="assembly") assembly(world,atoi(arg(argc,argv,"--N","173")),atoi(arg(argc,argv,"--tile","64")),atoi(arg(argc,argv,"--drop-chunk","0")),atoi(arg(argc,argv,"--L","28")));
        else if(mode=="factor" || mode=="residual") real_factor(world,argc,argv);
        else if(mode=="monolithic") monolithic(world,argc,argv);
        else throw std::invalid_argument("unknown mode");
        CHECK(cudaStreamSynchronize(0));
        unsigned long long high=0,reserve=0; size_t free=0,total=0;
        CHECK(cudaMemPoolGetAttribute(pool,cudaMemPoolAttrUsedMemHigh,&high));
        CHECK(cudaMemPoolGetAttribute(pool,cudaMemPoolAttrReservedMemCurrent,&reserve));
        CHECK(cudaMemGetInfo(&free,&total));
        struct rusage u; getrusage(RUSAGE_SELF,&u);
        std::cout << "MEMORY rank="<<rank<<" pool_peak="<<high<<" outside_pool_now="<<(total-free-reserve)
          <<" host_peak="<<(uint64_t(u.ru_maxrss)*1024)<<" wall="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<"\n";
      }
      MPI_Finalize(); return 0;
    }
  catch(const std::exception &e)
    { std::cerr<<"rank "<<rank<<": "<<e.what()<<"\n"; MPI_Abort(MPI_COMM_WORLD,1); return 1; }
}
