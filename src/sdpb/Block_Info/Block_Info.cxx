#include "../Block_Info.hxx"

Block_Info::Block_Info(const boost::filesystem::path &sdp_directory,
                       const boost::filesystem::path &checkpoint_in,
                       const size_t &procs_per_node,
                       const size_t &proc_granularity,
                       const Verbosity &verbosity)
{
  read_block_info(sdp_directory, checkpoint_in);
  std::vector<Block_Cost> block_costs(
    read_block_costs(sdp_directory, checkpoint_in));
  allocate_blocks(block_costs, procs_per_node, proc_granularity, verbosity);
}
