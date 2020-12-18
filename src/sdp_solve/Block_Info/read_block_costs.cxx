#include "../Block_Info.hxx"

#include <boost/filesystem/fstream.hpp>

std::vector<Block_Cost>
Block_Info::read_block_costs(const boost::filesystem::path &sdp_directory,
                             const boost::filesystem::path &checkpoint_in)
{
  const boost::filesystem::path sdp_block_timings_path(sdp_directory
                                                       / "block_timings"),
    checkpoint_block_timings_path(checkpoint_in / "block_timings");

  if(exists(checkpoint_in / ("checkpoint." + std::to_string(El::mpi::Rank()))))
    {
      if(exists(checkpoint_block_timings_path))
        {
          block_timings_filename = checkpoint_block_timings_path;
        }
    }
  else
    {
      block_timings_filename
        = (exists(checkpoint_block_timings_path)
             ? checkpoint_block_timings_path
             : (exists(sdp_block_timings_path) ? sdp_block_timings_path : ""));
    }

  std::vector<Block_Cost> result;
  if(!block_timings_filename.empty())
    {
      size_t index(0), cost;
      boost::filesystem::ifstream costs(block_timings_filename);
      costs >> cost;
      while(costs.good())
        {
          result.emplace_back(cost, index);
          ++index;
          costs >> cost;
        }
      if(result.size() != num_points.size())
        {
          throw std::runtime_error(
            "Incompatible number of entries in '"
            + block_timings_filename.string() + "'. Expected "
            + std::to_string(num_points.size()) + " but found "
            + std::to_string(result.size()));
        }
    }
  else
    {
      // If no information, assign a cost proportional to the matrix
      // size.  This should balance out memory use when doing a timing
      // run.
      auto schur_sizes(schur_block_sizes());
      for(size_t block = 0; block < schur_sizes.size(); ++block)
        {
          result.emplace_back(schur_sizes[block] * schur_sizes[block], block);
        }
    }
  return result;
}
