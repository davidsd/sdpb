#include "../Block_Info.hxx"

#include "../../compute_block_grid_mapping.hxx"

#include <boost/filesystem/fstream.hpp>

namespace
{
  void
  read_vector_with_index(std::ifstream &input_stream,
                         const std::vector<size_t> &indices,
                         const size_t &index_scale, std::vector<size_t> &v)
  {
    std::vector<size_t> file_v;
    read_vector(input_stream, file_v);
    for(size_t index = 0; index < file_v.size(); ++index)
      {
        const size_t mapped_index(index_scale * indices[index / index_scale]
                                  + index % index_scale);
        if(mapped_index >= v.size())
          {
            v.resize(mapped_index + 1);
          }
        v[mapped_index] = file_v[index];
      }
  }
}

void Block_Info::read_block_info(const boost::filesystem::path &sdp_directory)
{
  size_t file_rank(0);
  do
    {
      const boost::filesystem::path block_path(
        sdp_directory / ("blocks." + std::to_string(file_rank)));
      boost::filesystem::ifstream block_stream(block_path);
      if(!block_stream.good())
        {
          throw std::runtime_error("Could not open '" + block_path.string()
                                   + "'");
        }
      block_stream >> file_num_procs;
      if(!block_stream.good())
        {
          throw std::runtime_error("Corrupted file: " + block_path.string());
        }
      file_block_indices.emplace_back();
      auto &file_block_index(file_block_indices.back());
      read_vector(block_stream, file_block_index);

      read_vector_with_index(block_stream, file_block_index, 1, dimensions);
      read_vector_with_index(block_stream, file_block_index, 1, degrees);
      read_vector_with_index(block_stream, file_block_index, 1,
                             schur_block_sizes);
      read_vector_with_index(block_stream, file_block_index, 2,
                             psd_matrix_block_sizes);
      read_vector_with_index(block_stream, file_block_index, 2,
                             bilinear_pairing_block_sizes);
      ++file_rank;
    }
  while(file_rank < file_num_procs);
}
