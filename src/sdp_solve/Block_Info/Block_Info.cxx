#include "../Block_Info.hxx"

Block_Info::Block_Info(const boost::filesystem::path &sdp_directory,
                       const boost::filesystem::path &checkpoint_in,
                       const size_t &procs_per_node,
                       const size_t &proc_granularity,
                       const Verbosity &verbosity)
{
  read_block_info(sdp_directory);
  std::vector<Block_Cost> block_costs(
    read_block_costs(sdp_directory, checkpoint_in));
  allocate_blocks(block_costs, procs_per_node, proc_granularity, verbosity);
}

Block_Info::Block_Info(const boost::filesystem::path &sdp_directory,
                       const El::Matrix<int32_t> &block_timings,
                       const size_t &procs_per_node,
                       const size_t &proc_granularity,
                       const Verbosity &verbosity)
{
  read_block_info(sdp_directory);
  std::vector<Block_Cost> block_costs;
  for(int64_t block = 0; block < block_timings.Height(); ++block)
    {
      block_costs.emplace_back(block_timings(block, 0), block);
    }
  allocate_blocks(block_costs, procs_per_node, proc_granularity, verbosity);
}

Block_Info::Block_Info(const std::vector<size_t> &matrix_dimensions,
                       const size_t &procs_per_node,
                       const size_t &proc_granularity,
                       const Verbosity &verbosity)
    // TODO: This does not set the filename, file_block_indices, or
    // file_num_procs, since those are only useful when reading in info
    // from a filesystem.
    : dimensions(matrix_dimensions), degrees(matrix_dimensions.size(), 0),
      schur_block_sizes(matrix_dimensions.size()),
      psd_matrix_block_sizes(matrix_dimensions.size() * 2),
      bilinear_pairing_block_sizes(matrix_dimensions.size() * 2)
{
  for(size_t index(0); index != dimensions.size(); ++index)
    {
      schur_block_sizes.at(index)
        = dimensions.at(index) * (dimensions.at(index) + 1) / 2;
      psd_matrix_block_sizes.at(2 * index) = dimensions.at(index);
      psd_matrix_block_sizes.at(2 * index + 1) = 0;

      bilinear_pairing_block_sizes.at(2 * index)
        = bilinear_pairing_block_sizes.at(2 * index + 1)
        = dimensions.at(index);
    }

  std::vector<Block_Cost> block_costs;
  for(size_t block = 0; block < schur_block_sizes.size(); ++block)
    {
      block_costs.emplace_back(
        schur_block_sizes[block] * schur_block_sizes[block], block);
    }
  allocate_blocks(block_costs, procs_per_node, proc_granularity, verbosity);
}
