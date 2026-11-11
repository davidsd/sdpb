#include "Distributed_Q.hxx"
#include "Distributed_Q_Wide.cuh"
#include "../solver/Mpf_Ops.hxx"
#include "../rns/rns_kernels.cuh"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <sys/resource.h>

#define DQ_CHECK(x) RNS_CUDA_CHECK(x)

namespace sdpb_gpu
{
  namespace
  {
    constexpr int64_t zero_exp = -(int64_t(1) << 40);
    constexpr int threads = 256;
    int blocks(int64_t n) { return int((n + threads - 1) / threads); }
    int checked_dimension(int n,const char *name)
    {
      if(n<=0) throw std::invalid_argument(name);
      return n;
    }

    struct Tile
    {
      Mpf_Batch b;
      Tile() = default;
      Tile(int L, int rows, int cols)
      {
        b = {nullptr, nullptr, L, 1, rows, cols};
        DQ_CHECK(cudaMallocAsync(&b.m, size_t(L) * rows * cols * 4, 0));
        DQ_CHECK(cudaMallocAsync(&b.e, size_t(rows) * cols * 8, 0));
      }
      Tile(const Tile &) = delete;
      Tile &operator=(const Tile &) = delete;
      Tile(Tile &&o) noexcept : b(o.b) { o.b = {}; }
      Tile &operator=(Tile &&o) noexcept
      {
        if(this != &o) { clear(); b = o.b; o.b = {}; }
        return *this;
      }
      ~Tile() { clear(); }
      void clear()
      {
        if(b.m) cudaFreeAsync(b.m, 0);
        if(b.e) cudaFreeAsync(b.e, 0);
        b = {};
      }
    };
    template<class T> struct Buffer
    {
      T *p = nullptr;
      explicit Buffer(size_t n) { DQ_CHECK(cudaMallocAsync(&p, std::max<size_t>(n, 1) * sizeof(T), 0)); }
      ~Buffer() { if(p) cudaFreeAsync(p, 0); }
      Buffer(const Buffer &) = delete;
    };
    __global__ void copy_view(Mpf_Batch a, Mpf_Batch b)
    {
      const int64_t q=int64_t(blockIdx.x)*blockDim.x+threadIdx.x;
      const int64_t n=a.rows*a.cols;
      if(q>=n) return;
      const int64_t ai=(q/a.rows)*(a.ld?a.ld:a.rows)+q%a.rows;
      const int64_t bi=(q/b.rows)*(b.ld?b.ld:b.rows)+q%b.rows;
      const int64_t ap=a.plane?a.plane:n, bp=b.plane?b.plane:n;
      b.e[bi]=a.e[ai];
      for(int t=0;t<a.L;++t) b.m[int64_t(t)*bp+bi]=a.m[int64_t(t)*ap+ai];
    }
    void copy(const Mpf_Batch &a, Mpf_Batch &b)
    {
      if(a.rows != b.rows || a.cols != b.cols || a.L != b.L)
        throw std::invalid_argument("Distributed_Q copy shape");
      copy_view<<<blocks(a.entries()),threads>>>(a,b);
      DQ_CHECK(cudaGetLastError());
    }
    void broadcast(Mpf_Batch &a, Gpu_Collective &w, int root)
    {
      w.broadcast(a.m, size_t(a.L) * a.entries() * 4, root);
      w.broadcast(a.e, size_t(a.entries()) * 8, root);
    }
    void gemv(const Mpf_Batch &a, bool trans, const Mpf_Batch &v, Mpf_Batch &out)
    {
      Gemv_Spec s;
      s.no = trans ? a.cols : a.rows;
      s.nti = trans ? a.rows : a.cols;
      s.ao = trans ? a.lda() : 1;
      s.ati = trans ? 1 : a.lda();
      s.vti = 1;
      mpf_gemv(a, v, out, s, 0);
    }
    __global__ void diagonal_tops(const int64_t *e, uint64_t *tops, int b, int off)
    {
      int i = blockIdx.x * blockDim.x + threadIdx.x;
      if(i < b) tops[off+i] = uint64_t(e[int64_t(i)*b+i]) ^ (uint64_t(1)<<63);
    }
    __global__ void decode_scales(const uint64_t *tops, int64_t *s, int n, int L)
    {
      int i = blockIdx.x * blockDim.x + threadIdx.x;
      if(i < n) s[i] = (int64_t(tops[i] ^ (uint64_t(1)<<63)) + 32LL*L - 1) >> 1;
    }
    __global__ void scale_tile(int64_t *e, const int64_t *s, int b, int i0, int j0,
                               int mode)
    {
      const int64_t q = int64_t(blockIdx.x)*blockDim.x+threadIdx.x;
      if(q >= int64_t(b)*b || e[q] == zero_exp) return;
      const int i = i0+q%b, j = j0+q/b;
      // 0: Q -> D Q D; 1: Ls -> D^-1 Ls; 2: Xs -> Xs D
      e[q] += mode == 0 ? -s[i]-s[j] : (mode == 1 ? s[i] : -s[j]);
    }
    __global__ void pad_tile(int32_t *m, int64_t *e, int b, int L,
                             int i0, int j0, int N)
    {
      const int64_t q = int64_t(blockIdx.x)*blockDim.x+threadIdx.x;
      if(q >= int64_t(b)*b) return;
      const int i = i0+q%b, j = j0+q/b;
      if(i < N && j < N) return;
      const bool one = i == j;
      for(int t=0; t<L; ++t)
        m[int64_t(t)*b*b+q] = one && t == L-1 ? (1u<<29) : 0;
      e[q] = one ? -(32LL*L-3) : zero_exp;
    }
    __global__ void perturb(int64_t *e, int bit)
    {
      // A factor entry multiplied by two is an unambiguous negative control.
      // bit is a selected entry, not a requested error size.
      if(threadIdx.x == 0 && blockIdx.x == 0 && e[bit] != zero_exp) ++e[bit];
    }
    uint64_t digest(const void *data, size_t n)
    {
      auto p = static_cast<const unsigned char *>(data);
      uint64_t h=14695981039346656037ULL;
      for(size_t i=0;i<n;++i) { h ^= p[i]; h *= 1099511628211ULL; }
      return h;
    }

    struct Q_Pointer { const int32_t *m; const int64_t *e; };
    __device__ Q_Pointer q_entry(const Q_Pointer *tiles,int nt,int b,
                                 int row,int col,int64_t &entry)
    {
      int i=row/b,j=col/b;
      if(i>=j) { entry=int64_t(col%b)*b+row%b; return tiles[i*nt+j]; }
      entry=int64_t(row%b)*b+col%b;return tiles[j*nt+i];
    }
    __global__ void residual_top(const Q_Pointer *tiles,int nt,int b,int N,
                                  const int64_t *ve,const int64_t *rhs,
                                  uint64_t *top)
    {
      const int row=blockIdx.x*blockDim.x+threadIdx.x;
      if(row>=N)return;
      int64_t t=INT64_MIN;
      // Treat the L-word RHS as a product with an extra L zero low words.
      if(rhs[row]!=zero_exp)t=rhs[row]-32*dq_wide::L;
      for(int col=0;col<N;++col)
        {
          int64_t q;auto a=q_entry(tiles,nt,b,row,col,q);
          if(!a.e || a.e[q]==zero_exp || ve[col]==zero_exp)continue;
          const int64_t exponent=a.e[q]+ve[col];
          if(exponent>t)t=exponent;
        }
      top[row]=uint64_t(t)^(uint64_t(1)<<63);
    }
    __global__ void residual_parts(const Q_Pointer *tiles,int nt,int b,int N,
                                    Mpf_Batch v,const uint64_t *top,int slices,
                                    uint32_t *part)
    {
      using namespace dq_wide;
      const int64_t tid=int64_t(blockIdx.x)*blockDim.x+threadIdx.x;
      const int row=int(tid%N),sl=int(tid/N);
      if(sl>=slices)return;
      u32 acc[W]={};
      const i64 base=top[row]?i64(top[row]^(uint64_t(1)<<63))-guard:0;
      const i64 vp=v.plane?v.plane:N;
      if(top[row])for(int col=sl;col<N;col+=slices)
        {
          i64 q;auto a=q_entry(tiles,nt,b,row,col,q);
          if(!a.e || a.e[q]==zero_exp || v.e[col]==zero_exp)continue;
          u32 x[L],y[L],p[P];
          const bool nx=magnitude(a.m,int64_t(b)*b,q,x);
          const bool ny=magnitude(v.m,vp,col,y);
          product(x,y,p);
          add_shifted(acc,p,a.e[q]+v.e[col]-base,nx!=ny);
        }
      for(int k=0;k<W;++k)part[(int64_t(sl)*W+k)*N+row]=acc[k];
    }
    __global__ void residual_local_sum(const uint32_t *part,int N,int slices,uint32_t *sum)
    {
      using namespace dq_wide;
      const int row=blockIdx.x*blockDim.x+threadIdx.x;
      if(row>=N)return;
      u32 acc[W]={};
      for(int sl=0;sl<slices;++sl)add_words(acc,part+int64_t(sl)*W*N+row,N);
      for(int k=0;k<W;++k)sum[int64_t(k)*N+row]=acc[k];
    }
    __global__ void residual_finish(const uint32_t *gather,int G,int N,
                                     const uint64_t *top,Mpf_Batch rhs,Mpf_Batch out)
    {
      using namespace dq_wide;
      const int row=blockIdx.x*blockDim.x+threadIdx.x;
      if(row>=N)return;
      u32 acc[W]={},r[L];i64 ex=zero;
      if(top[row])
        {
          for(int rank=0;rank<G;++rank)add_words(acc,gather+int64_t(rank)*W*N+row,N);
          negate(acc); // rhs - Q*v, before any rounding or floating reduction
          const i64 base=i64(top[row]^(uint64_t(1)<<63))-guard;
          if(rhs.e[row]!=zero_exp)
            {
              u32 x[L];const bool negative=magnitude(rhs.m,rhs.plane?rhs.plane:N,row,x);
              add_shifted(acc,x,rhs.e[row]-base,negative);
            }
          rounded(acc,base,r,ex);
        }
      else for(int k=0;k<L;++k)r[k]=0;
      const i64 plane=out.plane?out.plane:N;
      for(int k=0;k<L;++k)out.m[int64_t(k)*plane+row]=int32_t(r[k]);
      out.e[row]=ex;
    }
  }

  void distributed_sum(const Mpf_Batch &local, Mpf_Batch &sum,
                        Gpu_Collective &w)
  {
    if(local.batch != 1 || local.ld || local.plane || local.bstride)
      throw std::invalid_argument("distributed_sum requires a packed matrix");
    if(w.size() == 1) { copy(local, sum); return; }
    const int64_t n = local.entries();
    Buffer<int32_t> m(size_t(w.size())*local.L*n);
    Buffer<int64_t> e(size_t(w.size())*n);
    w.allgather(local.m, m.p, size_t(local.L)*n*4);
    w.allgather(local.e, e.p, size_t(n)*8);
    mpf_zero(sum, 0);
    for(int r=0;r<w.size();++r)
      {
        Mpf_Batch a{m.p+size_t(r)*local.L*n, e.p+size_t(r)*n,
                    local.L, 1, local.rows, local.cols};
        mpf_batch_add(sum, a, 1, sum, 0);
      }
  }

  std::string full_Q_tile_digests(const Mpf_Batch &Q,int b)
  {
    if(Q.batch!=1 || Q.rows!=Q.cols || b<=0 || b%32)
      throw std::invalid_argument("full_Q_tile_digests dimensions");
    const int N=int(Q.rows),T=(N+b-1)/b;
    Tile tile(Q.L,b,b);
    std::vector<int32_t> m(size_t(Q.L)*b*b);
    std::vector<int64_t> e(size_t(b)*b);
    std::ostringstream out;
    for(int j=0;j<T;++j) for(int i=j;i<T;++i)
      {
        mpf_zero(tile.b,0);
        const int rows=std::min(b,N-i*b),cols=std::min(b,N-j*b);
        auto src=Q.block(i*b,j*b,rows,cols);
        auto dest=tile.b.block(0,0,rows,cols);
        copy(src,dest);
        pad_tile<<<blocks(int64_t(b)*b),threads>>>(tile.b.m,tile.b.e,b,Q.L,i*b,j*b,N);
        DQ_CHECK(cudaGetLastError());
        mpf_download(tile.b,m.data(),e.data(),0);
        out << i << " " << j << " " << std::hex << digest(m.data(),m.size()*4)
            << " " << digest(e.data(),e.size()*8) << std::dec << "\n";
      }
    return out.str();
  }

  struct Distributed_Q::Impl
  {
    int N, L, D, b, nt, np;
    int64_t rows;
    Gpu_Collective &world;
    Gram_Options go;
    Distributed_Q_Options opt;
    std::vector<std::pair<int,int>> order;
    std::vector<int> index;
    std::vector<Tile> Q, F, Xdiag;
    std::unique_ptr<Block_Factor> engine;
    std::map<std::string,double> seconds;
    uint64_t peak_device = 0, peak_host = 0, operation_peak = 0;
    using Clock = std::chrono::steady_clock;
    Impl(int n, int l, int d, int64_t r, Gpu_Collective &w,
         Gram_Options g, Distributed_Q_Options o)
      : N(checked_dimension(n,"Distributed_Q N")), L(l), D(d),
        b(checked_dimension(o.tile,"Distributed_Q tile")), nt((N+b-1)/b), np(nt*b), rows(r),
        world(w), go(g), opt(o), index(size_t(nt)*nt, -1), Xdiag(nt)
    {
      if(L!=28 || D<=0 || rows<=0 || b%32 || opt.row_panel<=0 || opt.row_panel%16
         || opt.row_panel > (1<<17))
        throw std::invalid_argument("Distributed_Q dimensions");
      // Packed tile columns, distributed cyclically by tile (not by column).
      for(int j=0;j<nt;++j) for(int i=j;i<nt;++i)
        {
          index[size_t(i)*nt+j] = int(order.size());
          order.emplace_back(i,j);
        }
      Q.resize(order.size()); F.resize(order.size());
      for(size_t t=0;t<order.size();++t)
        if(int(t%w.size()) == w.rank()) Q[t] = Tile(L,b,b);
      engine = std::make_unique<Block_Factor>(32,1,L,D,go,cudaStream_t(0));
    }
    int id(int i,int j) const
    {
      if(j<0 || i<j || i>=nt) throw std::out_of_range("Q tile index");
      return index[size_t(i)*nt+j];
    }
    int owner(int i,int j) const { return id(i,j)%world.size(); }
    bool local(int i,int j) const { return owner(i,j)==world.rank(); }
    void done(const char *name, Clock::time_point t)
    {
      DQ_CHECK(cudaStreamSynchronize(0));
      seconds[name] += std::chrono::duration<double>(Clock::now()-t).count();
      int device=0; cudaMemPool_t pool; size_t free=0,total=0;
      DQ_CHECK(cudaGetDevice(&device));
      DQ_CHECK(cudaDeviceGetDefaultMemPool(&pool,device));
      uint64_t high=0,reserved=0,zero=0;
      DQ_CHECK(cudaMemPoolGetAttribute(pool,cudaMemPoolAttrUsedMemHigh,&high));
      DQ_CHECK(cudaMemPoolGetAttribute(pool,cudaMemPoolAttrReservedMemCurrent,&reserved));
      DQ_CHECK(cudaMemGetInfo(&free,&total));
      const uint64_t observed=high+uint64_t(total-free>=reserved?total-free-reserved:0);
      peak_device=std::max(peak_device,observed);
      operation_peak=std::max(operation_peak,observed);
      DQ_CHECK(cudaMemPoolSetAttribute(pool,cudaMemPoolAttrUsedMemHigh,&zero));
      struct rusage u; getrusage(RUSAGE_SELF,&u);
      peak_host=std::max(peak_host,uint64_t(u.ru_maxrss)*1024);
    }
    void matmul(const Mpf_Batch &a,bool ta,const Mpf_Batch &c,bool tc,Mpf_Batch &out)
    { engine->mm(a,ta,c,tc,out); }
    Tile get(const std::vector<Tile> &tiles, int i,int j)
    {
      Tile t(L,b,b);
      if(local(i,j)) copy(tiles[id(i,j)].b,t.b);
      broadcast(t.b,world,owner(i,j));
      return t;
    }
    void triangular(const Mpf_Batch &rhs, Mpf_Batch &out, bool trans)
    {
      if(rhs.rows != N || rhs.cols != 1 || out.rows != N || out.cols != 1)
        throw std::invalid_argument("distributed Q solve vector shape");
      Tile state(L,np,1), result(L,np,1);
      mpf_zero(state.b,0); mpf_zero(result.b,0);
      auto in = state.b.block(0,0,N,1);
      copy(rhs,in);
      for(int step=0;step<nt;++step)
        {
          const int k=trans ? nt-1-step : step;
          Tile v(L,b,1), z(L,b,1), temp(L,b,1), residual(L,b,1), correction(L,b,1);
          auto sk=state.b.block(k*b,0,b,1);
          copy(sk,v.b);
          const int root=owner(k,k);
          if(world.rank()==root)
            {
              gemv(Xdiag[k].b,trans,v.b,z.b);
              // Refine diagonal substitution, using the actual triangular factor.
              gemv(F[id(k,k)].b,trans,z.b,temp.b);
              mpf_batch_add(v.b,temp.b,-1,residual.b,0);
              gemv(Xdiag[k].b,trans,residual.b,correction.b);
              mpf_batch_add(z.b,correction.b,1,z.b,0);
            }
          broadcast(z.b,world,root);
          auto rk=result.b.block(k*b,0,b,1);
          copy(z.b,rk);
          Tile partial(L,np,1), update(L,np,1);
          mpf_zero(partial.b,0);
          if(!trans)
            for(int i=k+1;i<nt;++i)
              if(local(i,k))
                {
                  gemv(F[id(i,k)].b,false,z.b,temp.b);
                  auto u=partial.b.block(i*b,0,b,1);
                  copy(temp.b,u);
                }
          if(trans)
            for(int j=0;j<k;++j)
              if(local(k,j))
                {
                  gemv(F[id(k,j)].b,true,z.b,temp.b);
                  auto u=partial.b.block(j*b,0,b,1);
                  copy(temp.b,u);
                }
          distributed_sum(partial.b,update.b,world);
          mpf_batch_add(state.b,update.b,-1,state.b,0);
        }
      auto cropped=result.b.block(0,0,N,1);
      copy(cropped,out);
    }
  };

  Distributed_Q::Distributed_Q(int N,int L,int digits,int64_t rows,
                               Gpu_Collective &w,Gram_Options go,
                               Distributed_Q_Options opt)
    : p_(new Impl(N,L,digits,rows,w,go,opt)) {}
  Distributed_Q::~Distributed_Q() = default;
  int Distributed_Q::tile_size() const { return p_->b; }
  int Distributed_Q::tiles() const { return p_->nt; }
  int Distributed_Q::owner(int i,int j) const { return p_->owner(i,j); }

  void Distributed_Q::assemble(const std::vector<Mpf_Matrix> &P)
  {
    auto &s=*p_; const int b=s.b, G=s.world.size();
    s.operation_peak=0;
    auto t=Impl::Clock::now();
    Buffer<uint64_t> tops(s.N);
    DQ_CHECK(cudaMemsetAsync(tops.p,0,size_t(s.N)*8,0));
    for(const auto &a:P)
      {
        if(a.cols != s.N || a.L != s.L)
          throw std::invalid_argument("Distributed_Q P shape");
        mpf_scan_column_tops(a,tops.p);
      }
    s.world.max_u64(tops.p,s.N);
    s.done("column maxima",t);
    Gram_Options gd=s.go;
    gd.max_rect=b; gd.base=std::min(256,b);
    gd.crt_batch=std::min<int64_t>(int64_t(b)*b,gd.crt_batch);
    gd.tile_rows=gd.tile_cols=0;
    Mpf_Gram diag(b,s.L,s.D,s.rows,gd);
    gd.tile_rows=gd.tile_cols=b;
    Mpf_Gram rect(2*b,s.L,s.D,s.rows,gd);
    const int r=int(rect.gram().basis().primes.size());
    if(diag.gram().basis().primes != rect.gram().basis().primes
       || uint64_t(G)*(rect.gram().basis().primes.front()-1)>UINT32_MAX)
      throw std::runtime_error("Distributed_Q residue collective overflow/basis");
    const size_t count=size_t(r)*b*b;
    Buffer<uint32_t> send(size_t(G)*count), recv(count);
    Buffer<int32_t> map(2*b);
    std::vector<int32_t> hm(2*b);
    auto configure=[&](int tile) -> Mpf_Gram & {
      auto [i,j]=s.order[tile];
      const bool diagonal=i==j;
      for(int x=0;x<b;++x)
        {
          hm[x]=j*b+x<s.N ? j*b+x : -1;
          hm[b+x]=i*b+x<s.N ? i*b+x : -1;
        }
      DQ_CHECK(cudaMemcpyAsync(map.p,hm.data(),size_t(diagonal?b:2*b)*4,cudaMemcpyHostToDevice,0));
      auto &g=diagonal?diag:rect;
      g.set_column_tops(tops.p,map.p);
      // The host map buffer must not be changed until this copy completes.
      DQ_CHECK(cudaStreamSynchronize(0));
      return g;
    };
    for(size_t start=0;start<s.order.size();start+=G)
      {
        DQ_CHECK(cudaMemsetAsync(send.p,0,size_t(G)*count*4,0));
        t=Impl::Clock::now();
        for(int dest=0;dest<G && start+dest<s.order.size();++dest)
          {
            auto &g=configure(int(start+dest));
            g.gram().reset();
            for(const auto &a:P)
              for(int64_t row=0;row<a.rows;row+=s.opt.row_panel)
                {
                  auto panel=a;
                  panel.m+=row; panel.e+=row;
                  panel.rows=std::min<int64_t>(s.opt.row_panel,a.rows-row);
                  panel.ld=a.ld?a.ld:a.rows;
                  panel.plane=a.plane?a.plane:a.rows*a.cols;
                  g.add_columns(panel,map.p);
                }
            g.gram().export_residues(send.p+size_t(dest)*count);
          }
        s.done("residues and products",t);
        t=Impl::Clock::now();
        s.world.reduce_scatter_u32(send.p,recv.p,count);
        s.done("residue reduce-scatter",t);
        const size_t my=start+s.world.rank();
        t=Impl::Clock::now();
        if(my<s.order.size())
          {
            auto &g=configure(int(my));
            g.gram().import_residues(recv.p);
            auto [i,j]=s.order[my];
            if(i==j)
              {
                Tile packed(s.L,int(int64_t(b)*(b+1)/2),1);
                g.finish(packed.b.m,packed.b.e);
                mpf_unpack_symmetric(packed.b,s.Q[my].b,0);
              }
            else g.finish(s.Q[my].b.m,s.Q[my].b.e);
            pad_tile<<<blocks(int64_t(b)*b),threads>>>(s.Q[my].b.m,s.Q[my].b.e,b,s.L,i*b,j*b,s.N);
            DQ_CHECK(cudaGetLastError());
          }
        s.done("owner CRT",t);
      }
  }

  void Distributed_Q::load_tile(int i,int j,const int32_t *m,const int64_t *e)
  {
    auto &s=*p_;
    if(!s.local(i,j)) throw std::invalid_argument("load_tile on nonowner");
    auto &a=s.Q[s.id(i,j)].b;
    mpf_upload(m,e,a,0);
    pad_tile<<<blocks(a.entries()),threads>>>(a.m,a.e,s.b,s.L,i*s.b,j*s.b,s.N);
    DQ_CHECK(cudaGetLastError());
  }

  void Distributed_Q::download_tile(int i,int j,int32_t *m,int64_t *e,bool factor)
  {
    auto &s=*p_;
    if(!s.local(i,j)) throw std::invalid_argument("download_tile on nonowner");
    mpf_download((factor?s.F:s.Q)[s.id(i,j)].b,m,e,0);
  }

  void Distributed_Q::factor()
  {
    auto &s=*p_; const int b=s.b;
    s.operation_peak=0;
    auto t=Impl::Clock::now();
    Buffer<uint64_t> tops(s.np); Buffer<int64_t> scale(s.np);
    DQ_CHECK(cudaMemsetAsync(tops.p,0,size_t(s.np)*8,0));
    for(int k=0;k<s.nt;++k) if(s.local(k,k))
      diagonal_tops<<<blocks(b),threads>>>(s.Q[s.id(k,k)].b.e,tops.p,b,k*b);
    s.world.max_u64(tops.p,s.np);
    decode_scales<<<blocks(s.np),threads>>>(tops.p,scale.p,s.np,s.L);
    DQ_CHECK(cudaGetLastError());
    for(size_t q=0;q<s.order.size();++q) if(int(q%s.world.size())==s.world.rank())
      {
        s.F[q]=Tile(s.L,b,b); copy(s.Q[q].b,s.F[q].b);
        auto [i,j]=s.order[q];
        scale_tile<<<blocks(int64_t(b)*b),threads>>>(s.F[q].b.e,scale.p,b,i*b,j*b,0);
      }
    s.done("factor scaling",t);
    for(int k=0;k<s.nt;++k)
      {
        t=Impl::Clock::now();
        Tile diagonal(s.L,b,b), inverse(s.L,b,b);
        const int root=s.owner(k,k);
        if(s.world.rank()==root)
          {
            Factor_Options fo; fo.panel_refine=1; fo.correction_digits=-1;
            s.engine->factor(s.F[s.id(k,k)].b,diagonal.b,inverse.b,fo);
            copy(diagonal.b,s.F[s.id(k,k)].b);
            s.Xdiag[k]=Tile(s.L,b,b); copy(inverse.b,s.Xdiag[k].b);
          }
        broadcast(diagonal.b,s.world,root); broadcast(inverse.b,s.world,root);
        s.done("factor diagonal and broadcast",t);
        t=Impl::Clock::now();
        std::vector<Tile> panels(s.nt);
        for(int i=k+1;i<s.nt;++i)
          {
            panels[i]=Tile(s.L,b,b);
            if(s.local(i,k))
              {
                auto &a=s.F[s.id(i,k)].b;
                Tile product(s.L,b,b), residual(s.L,b,b), correction(s.L,b,b);
                s.matmul(a,false,inverse.b,true,panels[i].b);
                s.matmul(panels[i].b,false,diagonal.b,true,product.b);
                mpf_batch_add(a,product.b,-1,residual.b,0);
                s.matmul(residual.b,false,inverse.b,true,correction.b);
                mpf_batch_add(panels[i].b,correction.b,1,panels[i].b,0);
                copy(panels[i].b,a);
              }
            broadcast(panels[i].b,s.world,s.owner(i,k));
          }
        s.done("factor panel and broadcast",t);
        t=Impl::Clock::now();
        Tile product(s.L,b,b);
        for(int j=k+1;j<s.nt;++j) for(int i=j;i<s.nt;++i)
          if(s.local(i,j))
            {
              s.matmul(panels[i].b,false,panels[j].b,true,product.b);
              auto &a=s.F[s.id(i,j)].b;
              mpf_batch_add(a,product.b,-1,a,0,i==j);
            }
        s.engine->check();
        s.done("factor trailing update",t);
      }
    t=Impl::Clock::now();
    for(size_t q=0;q<s.order.size();++q) if(int(q%s.world.size())==s.world.rank())
      {
        auto [i,j]=s.order[q];
        scale_tile<<<blocks(int64_t(b)*b),threads>>>(s.F[q].b.e,scale.p,b,i*b,j*b,1);
        if(i==j)
          scale_tile<<<blocks(int64_t(b)*b),threads>>>(s.Xdiag[i].b.e,scale.p,b,i*b,j*b,2);
      }
    DQ_CHECK(cudaGetLastError());
    s.done("factor unscale",t);
  }

  void Distributed_Q::multiply(const Mpf_Batch &v,Mpf_Batch &out)
  {
    auto &s=*p_; const int b=s.b;
    Tile padded(s.L,s.np,1), local(s.L,s.np,1), sum(s.L,s.np,1), temp(s.L,b,1);
    mpf_zero(padded.b,0); mpf_zero(local.b,0);
    auto pv=padded.b.block(0,0,s.N,1); copy(v,pv);
    Tile vi(s.L,b,1),vj(s.L,b,1);
    for(size_t q=0;q<s.order.size();++q) if(int(q%s.world.size())==s.world.rank())
      {
        auto [i,j]=s.order[q];
        auto svi=padded.b.block(i*b,0,b,1), svj=padded.b.block(j*b,0,b,1);
        copy(svi,vi.b); copy(svj,vj.b);
        auto yi=local.b.block(i*b,0,b,1), yj=local.b.block(j*b,0,b,1);
        gemv(s.Q[q].b,false,vj.b,temp.b);
        mpf_batch_add(yi,temp.b,1,yi,0);
        if(i!=j)
          {
            gemv(s.Q[q].b,true,vi.b,temp.b);
            mpf_batch_add(yj,temp.b,1,yj,0);
          }
      }
    distributed_sum(local.b,sum.b,s.world);
    auto cropped=sum.b.block(0,0,s.N,1); copy(cropped,out);
  }

  void Distributed_Q::forward(const Mpf_Batch &rhs,Mpf_Batch &out)
  { p_->triangular(rhs,out,false); }

  void Distributed_Q::residual(const Mpf_Batch &rhs,const Mpf_Batch &v,Mpf_Batch &out)
  {
    auto &s=*p_;
    for(const auto *a:{&rhs,&v,static_cast<const Mpf_Batch *>(&out)})
      if(a->L!=28 || a->batch!=1 || a->rows!=s.N || a->cols!=1)
        throw std::invalid_argument("distributed Q residual vector shape");
    constexpr int slices=16,W=dq_wide::W;
    std::vector<Q_Pointer> host(size_t(s.nt)*s.nt,{nullptr,nullptr});
    for(size_t q=0;q<s.order.size();++q)if(int(q%s.world.size())==s.world.rank())
      { auto [i,j]=s.order[q];host[i*s.nt+j]={s.Q[q].b.m,s.Q[q].b.e}; }
    Buffer<Q_Pointer> tiles(host.size());
    DQ_CHECK(cudaMemcpy(tiles.p,host.data(),host.size()*sizeof(Q_Pointer),cudaMemcpyHostToDevice));
    Buffer<uint64_t> top(s.N);
    Buffer<uint32_t> part(size_t(slices)*W*s.N),local(size_t(W)*s.N),gather(size_t(s.world.size())*W*s.N);
    residual_top<<<blocks(s.N),threads>>>(tiles.p,s.nt,s.b,s.N,v.e,rhs.e,top.p);
    DQ_CHECK(cudaGetLastError());
    s.world.max_u64(top.p,s.N);
    residual_parts<<<blocks(int64_t(s.N)*slices),threads>>>(tiles.p,s.nt,s.b,s.N,v,top.p,slices,part.p);
    residual_local_sum<<<blocks(s.N),threads>>>(part.p,s.N,slices,local.p);
    DQ_CHECK(cudaGetLastError());
    s.world.allgather(local.p,gather.p,size_t(W)*s.N*4);
    residual_finish<<<blocks(s.N),threads>>>(gather.p,s.world.size(),s.N,top.p,rhs,out);
    DQ_CHECK(cudaGetLastError());
  }

  void Distributed_Q::solve(const Mpf_Batch &rhs,Mpf_Batch &out,int refine)
  {
    auto &s=*p_;
    s.operation_peak=0;
    if(refine<0) refine=s.opt.refine;
    auto t=Impl::Clock::now();
    Tile y(s.L,s.N,1);
    s.triangular(rhs,y.b,false); s.triangular(y.b,out,true);
    for(int step=0;step<refine;++step)
      {
        Tile ax(s.L,s.N,1), residual(s.L,s.N,1), correction(s.L,s.N,1);
        if(s.opt.wide_residual) this->residual(rhs,out,residual.b);
        else { multiply(out,ax.b); mpf_batch_add(rhs,ax.b,-1,residual.b,0); }
        s.triangular(residual.b,y.b,false); s.triangular(y.b,correction.b,true);
        mpf_batch_add(out,correction.b,1,out,0);
      }
    s.done("solve including refinement",t);
  }

  void Distributed_Q::factor_diagonal(Mpf_Batch &out)
  {
    auto &s=*p_; Tile local(s.L,s.np,1), sum(s.L,s.np,1), diag(s.L,s.b,1);
    mpf_zero(local.b,0);
    for(int i=0;i<s.nt;++i) if(s.local(i,i))
      {
        mpf_diagonal(s.F[s.id(i,i)].b,diag.b,0);
        auto dest=local.b.block(i*s.b,0,s.b,1);
        copy(diag.b,dest);
      }
    distributed_sum(local.b,sum.b,s.world);
    auto cropped=sum.b.block(0,0,s.N,1); copy(cropped,out);
  }

  std::string Distributed_Q::digests(bool factor)
  {
    auto &s=*p_; std::ostringstream out;
    std::vector<int32_t> m(size_t(s.L)*s.b*s.b);
    std::vector<int64_t> e(size_t(s.b)*s.b);
    for(size_t q=0;q<s.order.size();++q) if(int(q%s.world.size())==s.world.rank())
      {
        auto [i,j]=s.order[q];
        auto &a=(factor?s.F:s.Q)[q].b;
        mpf_download(a,m.data(),e.data(),0);
        out << i << " " << j << " " << std::hex << digest(m.data(),m.size()*4)
            << " " << digest(e.data(),e.size()*8) << std::dec << "\n";
      }
    return out.str();
  }

  void Distributed_Q::save_factor(const std::string &prefix)
  {
    auto &s=*p_;
    std::vector<int32_t> m(size_t(s.L)*s.b*s.b);
    std::vector<int64_t> e(size_t(s.b)*s.b);
    for(size_t q=0;q<s.order.size();++q) if(int(q%s.world.size())==s.world.rank())
      {
        auto [i,j]=s.order[q];
        mpf_download(s.F[q].b,m.data(),e.data(),0);
        auto name=prefix+"-"+std::to_string(i)+"-"+std::to_string(j);
        std::ofstream fm(name+".m.bin",std::ios::binary), fe(name+".e.bin",std::ios::binary);
        fm.write(reinterpret_cast<char *>(m.data()),m.size()*4);
        fe.write(reinterpret_cast<char *>(e.data()),e.size()*8);
        if(!fm || !fe) throw std::runtime_error("factor tile output failed");
      }
  }
  void Distributed_Q::corrupt_tile(int i,int j,int bit)
  {
    auto &s=*p_;
    if(bit<0 || bit>=s.b*s.b) throw std::invalid_argument("corrupt_tile entry");
    if(s.local(i,j)) perturb<<<1,1>>>(s.F[s.id(i,j)].b.e,bit);
    DQ_CHECK(cudaGetLastError());
  }
  uint64_t Distributed_Q::stored_bytes() const
  {
    auto &s=*p_; uint64_t n=0;
    for(auto *v:{&s.Q,&s.F,&s.Xdiag}) for(const auto &a:*v)
      if(a.b.m) n+=a.b.entries()*(4*s.L+8);
    return n;
  }
  uint64_t Distributed_Q::peak_device_bytes() const { return p_->peak_device; }
  uint64_t Distributed_Q::last_operation_peak_device_bytes() const { return p_->operation_peak; }
  std::string Distributed_Q::profile_json() const
  {
    auto &s=*p_; std::ostringstream out;
    out << "{\"rank\":" << s.world.rank() << ",\"G\":" << s.world.size()
        << ",\"N\":" << s.N << ",\"tile\":" << s.b
        << ",\"stored_bytes\":" << stored_bytes()
        << ",\"peak_device_bytes\":" << s.peak_device
        << ",\"peak_host_bytes\":" << s.peak_host
        << ",\"collective_received_payload\":" << s.world.received_bytes()
        << ",\"collective_calls\":" << s.world.collective_calls() << ",\"seconds\":{";
    bool first=true;
    for(auto &kv:s.seconds)
      { if(!first) out << ","; first=false; out << "\"" << kv.first << "\":" << std::setprecision(9) << kv.second; }
    out << "}}"; return out.str();
  }
}
