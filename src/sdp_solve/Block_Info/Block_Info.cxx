#include "../Block_Info.hxx"

namespace fs = std::filesystem;

Block_Info::Block_Info(const Environment &env, const fs::path &sdp_path,
                       const fs::path &checkpoint_in,
                       const size_t &proc_granularity,
                       const Verbosity &verbosity)
{
  read_block_info(sdp_path);
  std::vector<Block_Cost> block_costs(
    read_block_costs(sdp_path, checkpoint_in));
  allocate_blocks(env, block_costs, proc_granularity, verbosity);
}

Block_Info::Block_Info(const Environment &env, const fs::path &sdp_path,
                       const El::Matrix<int32_t> &block_timings,
                       const size_t &proc_granularity,
                       const Verbosity &verbosity)
{
  read_block_info(sdp_path);
  std::vector<Block_Cost> block_costs;
  for(int64_t block = 0; block < block_timings.Height(); ++block)
    {
      block_costs.emplace_back(block_timings(block, 0), block);
    }
  allocate_blocks(env, block_costs, proc_granularity, verbosity);
}

Block_Info::Block_Info(const Environment &env,
                       const std::vector<size_t> &matrix_dimensions,
                       const size_t &proc_granularity,
                       const Verbosity &verbosity)
    // TODO: This does not set the filename, file_block_indices, or
    // file_num_procs, since those are only useful when reading in info
    // from a filesystem.
    : dimensions(matrix_dimensions), num_points(matrix_dimensions.size(), 1)
{
  std::vector<Block_Cost> block_costs;
  auto schur_sizes(schur_block_sizes());
  for(size_t block = 0; block < schur_sizes.size(); ++block)
    {
      block_costs.emplace_back(schur_sizes[block] * schur_sizes[block], block);
    }
  allocate_blocks(env, block_costs, proc_granularity, verbosity);
}
