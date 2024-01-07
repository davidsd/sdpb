#include "PVM_SDP.hxx"

#include "PVM_File_Parse_Result.hxx"
#include "sdp_read/sdp_read.hxx"

namespace fs = std::filesystem;

namespace
{
  void validate(const PVM_SDP &sdp)
  {
    if(sdp.objective.empty())
      El::RuntimeError("Objectives not found in input files.");
  }
}

PVM_SDP::PVM_SDP(const fs::path &input_path) : PVM_SDP(std::vector{input_path})
{}

PVM_SDP::PVM_SDP(const std::vector<fs::path> &input_paths)
{
  num_matrices = 0;

  for(const auto &file : collect_files_expanding_nsv(input_paths))
    {
      // Simple round-robin for matrices across all files
      // TODO optimize https://github.com/davidsd/sdpb/issues/150
      auto should_parse_matrix = [this](size_t matrix_index) {
        return (this->num_matrices + matrix_index) % El::mpi::Size()
               == El::mpi::Rank();
      };
      PVM_File_Parse_Result result(file, should_parse_matrix);

      // TODO throw error if objectives are defined in several files?
      if(!result.objective.empty())
        objective = std::move(result.objective);

      for(auto &[index_in_file, matrix] : result.parsed_matrices)
        {
          size_t global_index = num_matrices + index_in_file;
          matrices.emplace_back(std::move(matrix));
          matrix_index_local_to_global.push_back(global_index);
        }

      num_matrices += result.num_matrices;
    }

  validate(*this);
}