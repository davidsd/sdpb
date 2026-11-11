#include "Block_Info.hxx"

#include "parse_sdpa/parse_sdpa.hxx"
#include "sdp_solve/Block_Info/block_timings.hxx"

#include <ostream>

namespace fs = std::filesystem;

namespace
{
  [[nodiscard]] std::vector<Block_Cost>
  block_dimensions_to_costs(const std::vector<size_t> &block_dimensions)
  {
    std::vector<Block_Cost> block_costs;
    block_costs.reserve(block_dimensions.size());
    for(size_t block_index = 0; block_index < block_dimensions.size();
        ++block_index)
      {
        const auto cost = block_dimensions.at(block_index)
                          * block_dimensions.at(block_index);
        block_costs.emplace_back(cost, block_index);
      }
    return block_costs;
  }

  [[nodiscard]] std::vector<Block_Cost>
  block_timings_to_costs(const El::Matrix<int32_t> &block_timings)
  {
    std::vector<Block_Cost> block_costs;
    for(int64_t block = 0; block < block_timings.Height(); ++block)
      {
        block_costs.emplace_back(block_timings(block, 0), block);
      }
    return block_costs;
  }

  [[nodiscard]] std::vector<Block_Cost> block_costs_from_timings_or_dimensions(
    const fs::path &block_timings_path,
    const std::vector<size_t> &block_dimensions)
  {
    if(!fs::exists(block_timings_path))
      return block_dimensions_to_costs(block_dimensions);
    return read_block_costs_from_timings(block_timings_path,
                                         block_dimensions.size());
  }
}

namespace Sdpb::Sdpa
{
  Block_Info::Block_Info(const Environment &env,
                         const size_t &primal_dimension,
                         const std::vector<size_t> &block_dimensions,
                         const std::vector<Block_Cost> &block_costs,
                         const size_t &proc_granularity,
                         const Verbosity &verbosity)
      : Abstract_Block_Info(
          env, block_costs, proc_granularity,
          [&block_dimensions](std::ostream &os,
                              const size_t block_index) -> auto & {
            return os << block_index << "(" << block_dimensions[block_index]
                      << ")";
          },
          verbosity),
        primal_dimension(primal_dimension),
        block_dimensions(block_dimensions)
  {
    ASSERT_EQUAL(block_dimensions.size(), block_costs.size());
  }

  Block_Info
  Block_Info::create(const Environment &env, const fs::path &sdp_path,
                     const fs::path &block_timings_path,
                     const size_t &proc_granularity,
                     const Verbosity &verbosity)
  {
    const auto structure = read_block_structure(sdp_path);
    const auto primal_dimension = structure.m_dim;
    const auto block_dimensions = structure.solver_block_sizes;
    const auto block_costs = block_costs_from_timings_or_dimensions(
      block_timings_path, block_dimensions);
    return Block_Info(env, primal_dimension, block_dimensions, block_costs,
                      proc_granularity, verbosity);
  }

  Block_Info
  Block_Info::create(const Environment &env, const fs::path &sdp_path,
                     const El::Matrix<int32_t> &block_timings,
                     const size_t &proc_granularity,
                     const Verbosity &verbosity)
  {
    const auto structure = read_block_structure(sdp_path);
    const auto primal_dimension = structure.m_dim;
    const auto block_dimensions = structure.solver_block_sizes;
    const auto block_costs = block_timings_to_costs(block_timings);
    return Block_Info(env, primal_dimension, block_dimensions, block_costs,
                      proc_granularity, verbosity);
  }

  Block_Info
  Block_Info::create(const Environment &env, const size_t &primal_dimension,
                     const std::vector<size_t> &block_dimensions,
                     const size_t &proc_granularity,
                     const Verbosity &verbosity)
  {
    const auto block_costs = block_dimensions_to_costs(block_dimensions);
    return Block_Info(env, primal_dimension, block_dimensions, block_costs,
                      proc_granularity, verbosity);
  }

  void swap(Block_Info &a, Block_Info &b) noexcept
  {
    using std::swap;
    swap(static_cast<Abstract_Block_Info &>(a),
         static_cast<Abstract_Block_Info &>(b));
    swap(a.primal_dimension, b.primal_dimension);
    swap(a.block_dimensions, b.block_dimensions);
  }
}
