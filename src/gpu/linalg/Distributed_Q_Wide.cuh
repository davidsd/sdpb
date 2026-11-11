// Internal exact-word arithmetic for distributed Q residuals.
// Matrix/vector storage remains L=28. A full product is shifted only AFTER
// multiplication; 128 low guard bits and 64 high headroom bits are retained.
// The global integer sum is rounded once, after subtracting it from the RHS.
#pragma once
#include <cstdint>
#ifdef __CUDACC__
#define DQW_HD __host__ __device__ __forceinline__
#define DQW_UNROLL _Pragma("unroll")
#else
#define DQW_HD inline
#define DQW_UNROLL
#endif
namespace sdpb_gpu { namespace dq_wide {
using u32=uint32_t;
using u64=uint64_t;
using i64=int64_t;
constexpr int L=28, P=2*L, guard=128, W=P+guard/32+2;
constexpr i64 zero=-(i64(1)<<40);

template<int K> DQW_HD void negate(u32 (&x)[K])
{
  u64 c=1;
  DQW_UNROLL
  for(int k=0;k<K;++k) { const u64 a=u64(~x[k])+c; x[k]=u32(a); c=a>>32; }
}
template<int K> DQW_HD bool magnitude(const int32_t *m,i64 plane,i64 index,u32 (&x)[K])
{
  DQW_UNROLL
  for(int k=0;k<K;++k) x[k]=u32(m[i64(k)*plane+index]);
  const bool negative=int32_t(x[K-1])<0;
  if(negative) negate(x);
  return negative;
}
DQW_HD void product(const u32 (&x)[L],const u32 (&y)[L],u32 (&p)[P])
{
  u32 c0=0,c1=0,c2=0;
  DQW_UNROLL
  for(int k=0;k<P-1;++k)
    {
      DQW_UNROLL
      for(int i=(k<L?0:k-L+1);i<L && i<=k;++i)
        {
#ifdef __CUDA_ARCH__
          asm("mad.lo.cc.u32 %0, %3, %4, %0;\n\t"
              "madc.hi.cc.u32 %1, %3, %4, %1;\n\t"
              "addc.u32 %2, %2, 0;"
              : "+r"(c0), "+r"(c1), "+r"(c2) : "r"(x[i]), "r"(y[k-i]));
#else
          const u64 v=u64(x[i])*y[k-i],lo=u64(c0)+u32(v);
          const u64 hi=u64(c1)+(v>>32)+(lo>>32);
          c0=u32(lo);c1=u32(hi);c2+=u32(hi>>32);
#endif
        }
      p[k]=c0;c0=c1;c1=c2;c2=0;
    }
  p[P-1]=c0;
}
// Word k of x*2^shift, with discarded low bits truncated toward zero.
template<int K> DQW_HD u32 shifted_word(const u32 (&x)[K],int k,i64 shift)
{
  const i64 bit=32LL*k-shift;
  const i64 q=bit>=0?bit/32:-((-bit+31)/32);
  const int r=int(bit-32*q);
  const u32 a=q>=0 && q<K?x[q]:0;
  const u32 b=q+1>=0 && q+1<K?x[q+1]:0;
  return (a>>r)|(r?b<<(32-r):0);
}
template<int K> DQW_HD void add_shifted(u32 (&acc)[W],const u32 (&x)[K],i64 shift,bool negative)
{
  const u32 mask=negative?~u32(0):0;
  u64 carry=negative?1:0;
  DQW_UNROLL
  for(int k=0;k<W;++k)
    {
      const u64 t=u64(acc[k])+(shifted_word(x,k,shift)^mask)+carry;
      acc[k]=u32(t);carry=t>>32;
    }
}
DQW_HD void add_words(u32 (&acc)[W],const u32 *p,i64 stride)
{
  u64 carry=0;
  DQW_UNROLL
  for(int k=0;k<W;++k)
    { const u64 t=u64(acc[k])+p[i64(k)*stride]+carry;acc[k]=u32(t);carry=t>>32; }
}
DQW_HD int word_bits(u32 x)
{
#ifdef __CUDA_ARCH__
  return x?32-__clz(x):0;
#else
  return x?32-__builtin_clz(x):0;
#endif
}
// Same normalized 894-bit magnitude and half-up rule as Mpf_Ops.
// The input is a signed W-word integer; result is a signed L-word integer.
DQW_HD void rounded(u32 (&acc)[W],i64 base,u32 (&out)[L],i64 &exponent)
{
  const bool negative=int32_t(acc[W-1])<0;
  if(negative) negate(acc);
  int bits=0;
  DQW_UNROLL
  for(int k=0;k<W;++k) if(acc[k]) bits=32*k+word_bits(acc[k]);
  if(!bits)
    { for(int k=0;k<L;++k) out[k]=0;exponent=zero;return; }
  const int shift=bits-(32*L-2);
  u64 carry=shift>0?(acc[(shift-1)/32]>>((shift-1)%32))&1:0;
  DQW_UNROLL
  for(int k=0;k<L;++k)
    { const u64 t=u64(shifted_word(acc,k,-i64(shift)))+carry;out[k]=u32(t);carry=t>>32; }
  exponent=base+shift;
  if(out[L-1]&0xc0000000u)
    {
      for(int k=0;k<L;++k) out[k]=0;
      out[L-1]=0x20000000u;++exponent;
    }
  if(negative) negate(out);
}
} }
#undef DQW_HD
#undef DQW_UNROLL
