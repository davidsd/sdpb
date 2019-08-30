#include "../../../read_vector.hxx"
#include "../../Block_Vector.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void read_primal_objective_c(const boost::filesystem::path &sdp_directory,
                             const std::vector<size_t> &block_indices,
                             const El::Grid &grid,
                             Block_Vector &primal_objective_c)
{
  primal_objective_c.blocks.reserve(block_indices.size());
  for(auto &block_index : block_indices)
    {
      const boost::filesystem::path primal_path(
        sdp_directory / ("primal_objective_c." + std::to_string(block_index)));
      boost::filesystem::ifstream primal_stream(primal_path);
      if(!primal_stream.good())
        {
          throw std::runtime_error("Could not open '" + primal_path.string()
                                   + "'");
        }

      std::vector<El::BigFloat> temp;
      read_vector(primal_stream, temp);
      primal_objective_c.blocks.emplace_back(temp.size(), 1, grid);
      auto &block(primal_objective_c.blocks.back());
      size_t local_height(block.LocalHeight());
      if(block.GlobalCol(0) == 0)
        {
          for(size_t row = 0; row < local_height; ++row)
            {
              size_t global_row(block.GlobalRow(row));
              block.SetLocal(row, 0, temp.at(global_row));
            }
        }
    }
}
