#include "../Block_Info.hxx"
#include "sdpb_util/Timers/Timers.hxx"

namespace fs = std::filesystem;

void read_objectives(const fs::path &sdp_path, const El::Grid &grid,
                     El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b,
                     Timers &timers);

std::vector<Block_Cost>
Block_Info::read_block_costs(const fs::path &sdp_path,
                             const fs::path &checkpoint_in,
                             const Environment &env)
{
  const fs::path sdp_block_timings_path(sdp_path / "block_timings"),
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
      std::ifstream costs(block_timings_filename);
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

      El::Grid grid(this->mpi_comm.value);
      El::BigFloat objective_const;
      El::DistMatrix<El::BigFloat> dual_objective_b;
      // TODO pass timers as argument
      Timers timers(env, false);
      // TODO objectives are already read in SDP::SDP(),
      // we should reuse them instead of reading again
      read_objectives(sdp_path, grid, objective_const, dual_objective_b,
                      timers);

      size_t dual_objective_size = dual_objective_b.Height(); // N
      auto schur_sizes = schur_block_sizes();
      auto psd_sizes = psd_matrix_block_sizes();
      auto bilinear_sizes = bilinear_pairing_block_sizes();

      auto elements_count
        = [](const std::vector<size_t> &sizes, const size_t index) {
            return sizes[index] * sizes[index];
          };

      for(size_t block = 0; block < schur_sizes.size(); ++block)
        {
          auto schur = elements_count(schur_sizes, block);
          auto psd = elements_count(psd_sizes, 2 * block)
                     + elements_count(psd_sizes, 2 * block + 1);
          auto bilinear = elements_count(bilinear_sizes, 2 * block)
                          + elements_count(bilinear_sizes, 2 * block + 1);
          // P'xN, a band of B(=free_var_matrix)
          auto B_band = schur_sizes[block] * dual_objective_size;
          // Estimate total RAM associated with the block.
          // (There is also a RAM contribution from #(Q)=NxN, but it's
          // block-independent)
          auto total_cost = 2 * B_band + 5 * psd + 2 * schur + 2 * bilinear;
          result.emplace_back(total_cost, block);
        }
    }
  return result;
}
