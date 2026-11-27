// Deterministic GPU chunks and whole-chunk ownership. No CUDA or MPI types.
// The chunk construction is the existing G=1 layout, including its ordering.
#pragma once
#include "Gpu_Solver.hxx"
#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

namespace sdpb_gpu
{
  inline std::vector<Chunk_Info>
  gpu_chunks(const std::vector<Block_Shape> &shapes, int cap, int64_t row_cap=0)
  {
    std::map<int, std::vector<int>> bydim;
    for(int j=0;j<int(shapes.size());++j) bydim[shapes[j].dim].push_back(j);
    std::vector<Chunk_Info> out;
    for(auto &kv:bydim)
      {
        auto &ids=kv.second;
        std::sort(ids.begin(),ids.end(),[&](int a,int b) {
          const auto &sa=shapes[a], &sb=shapes[b];
          const auto ka=std::make_pair(sa.K,std::max(sa.he,sa.ho));
          const auto kb=std::make_pair(sb.K,std::max(sb.he,sb.ho));
          return ka!=kb ? ka<kb : a<b;
        });
        const size_t step=cap>0 ? size_t(cap) : ids.size();
        for(size_t begin=0;begin<ids.size();)
          {
            Chunk_Info c{}; c.dim=kv.first;
            size_t i=begin;
            for(;i<std::min(ids.size(),begin+step);++i)
              {
                // A row cap ends a chunk early, at the block that would take it past
                // the cap; 0 keeps the blocks-only layout. Blocks are sorted by K, so
                // only the largest chunks of a dimension are cut.
                const int K=std::max(c.K,shapes[ids[i]].K);
                if(row_cap>0 && i>begin
                   && int64_t(i-begin+1)*(c.dim*(c.dim+1)/2*K)>row_cap) break;
                c.blocks.push_back(ids[i]);
                c.K=K;
                c.h=std::max({c.h,shapes[ids[i]].he,shapes[ids[i]].ho});
              }
            c.nb=int(c.blocks.size()); c.P=c.dim*(c.dim+1)/2*c.K;
            c.n=c.dim*c.h; out.push_back(std::move(c));
            begin=i;
          }
      }
    return out;
  }
  // Group blocks by padded matrix size. The row cap counts entries of
  // the padded matrices in one chunk of transformed constraints.
  inline std::vector<Chunk_Info>
  gpu_sdpa_chunks(const std::vector<int> &dims, int cap, int64_t row_cap)
  {
    constexpr int64_t round = 8;
    std::map<int, std::vector<int>> by_size;
    for(size_t j = 0; j < dims.size(); ++j)
      {
        const int64_t n = (int64_t(dims[j]) + round - 1) / round * round;
        if(dims[j] <= 0 || n * n > std::numeric_limits<int32_t>::max())
          throw std::invalid_argument("SDPA block exceeds the GPU index range");
        by_size[int(n)].push_back(int(j));
      }
    std::vector<Chunk_Info> out;
    for(const auto &[n, blocks] : by_size)
      {
        const int entries = n * n;
        // Matrix positions reserve their high bit for off-diagonal terms.
        // Transposing a batch uses one CUDA grid depth per block.
        size_t step = std::min(65535, std::numeric_limits<int32_t>::max() / entries);
        if(cap > 0)
          step = std::min(step, size_t(cap));
        if(row_cap > 0)
          step = std::min(step, size_t(std::max<int64_t>(1, row_cap / entries)));
        for(size_t begin = 0; begin < blocks.size(); begin += step)
          {
            Chunk_Info c{};
            c.dim = c.K = 1;
            c.h = c.n = n;
            c.P = entries;
            const size_t end = std::min(begin + step, blocks.size());
            c.blocks.assign(blocks.begin() + begin, blocks.begin() + end);
            c.nb = int(c.blocks.size());
            out.push_back(std::move(c));
          }
      }
    return out;
  }
  inline std::vector<int> gpu_chunk_owners(const std::vector<Chunk_Info> &chunks,int G)
  {
    if(G<1) throw std::invalid_argument("GPU rank count must be positive");
    std::vector<int> order(chunks.size()), owner(chunks.size());
    std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b) {
      return int64_t(chunks[a].nb)*chunks[a].P > int64_t(chunks[b].nb)*chunks[b].P;
    });
    std::vector<int64_t> load(G,0);
    for(int c:order)
      {
        const int rank=int(std::min_element(load.begin(),load.end())-load.begin());
        owner[c]=rank; load[rank]+=int64_t(chunks[c].nb)*chunks[c].P;
      }
    return owner;
  }
  inline std::vector<int> gpu_block_owners(const std::vector<Chunk_Info> &chunks,
                                           size_t num_blocks, int ranks)
  {
    const auto owners = gpu_chunk_owners(chunks, ranks);
    std::vector<int> block_owner(num_blocks, -1);
    for(size_t c = 0; c < chunks.size(); ++c)
      for(const int block : chunks[c].blocks)
        block_owner.at(block) = owners[c];
    return block_owner;
  }
}
