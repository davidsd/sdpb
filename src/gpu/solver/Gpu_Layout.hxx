// Deterministic GPU chunks and whole-chunk ownership. No CUDA or MPI types.
// The chunk construction is the existing G=1 layout, including its ordering.
#pragma once
#include "Gpu_Solver.hxx"
#include <algorithm>
#include <map>
#include <numeric>
#include <stdexcept>

namespace sdpb_gpu
{
  inline std::vector<Chunk_Info>
  gpu_chunks(const std::vector<Block_Shape> &shapes, int cap)
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
        for(size_t begin=0;begin<ids.size();begin+=step)
          {
            Chunk_Info c{}; c.dim=kv.first;
            for(size_t i=begin;i<std::min(ids.size(),begin+step);++i)
              {
                c.blocks.push_back(ids[i]);
                c.K=std::max(c.K,shapes[ids[i]].K);
                c.h=std::max({c.h,shapes[ids[i]].he,shapes[ids[i]].ho});
              }
            c.nb=int(c.blocks.size()); c.P=c.dim*(c.dim+1)/2*c.K;
            c.n=c.dim*c.h; out.push_back(std::move(c));
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
}
