#pragma once

#include "sdpb_util/block_mapping/Abstract_Block_Info.hxx"

#include <filesystem>

namespace Sdpb::Sdpa
{
  class Block_Info final : public Abstract_Block_Info
  {
  public:
    size_t primal_dimension;
    // Dimensions for all SDP blocks.
    std::vector<size_t> block_dimensions;

    [[nodiscard]]
    size_t num_blocks() const
    {
      return block_dimensions.size();
    }
    [[nodiscard]]
    size_t num_blocks_local() const
    {
      return block_indices.size();
    }

    // Helper functions to create Block_Info

    static Block_Info
    create(const Environment &env, const std::filesystem::path &sdp_path,
           const std::filesystem::path &block_timings_path,
           const size_t &proc_granularity, const Verbosity &verbosity);
    static Block_Info
    create(const Environment &env, const std::filesystem::path &sdp_path,
           const El::Matrix<int32_t> &block_timings,
           const size_t &proc_granularity, const Verbosity &verbosity);
    static Block_Info
    create(const Environment &env, const size_t &primal_dimension,
           const std::vector<size_t> &block_dimensions,
           const size_t &proc_granularity, const Verbosity &verbosity);

    Block_Info() = delete;

  private:
    Block_Info(const Environment &env, const size_t &primal_dimension,
               const std::vector<size_t> &block_dimensions,
               const std::vector<Block_Cost> &block_costs,
               const size_t &proc_granularity, const Verbosity &verbosity);
  };

  void swap(Block_Info &a, Block_Info &b) noexcept;
}
