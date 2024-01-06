#include "PMWP_SDP.hxx"

#include "PMWP_File_Parse_Result.hxx"
#include "sdp_read/sdp_read.hxx"

namespace fs = std::filesystem;

namespace
{
  void validate(const PMWP_SDP &sdp)
  {
    if(sdp.objective.empty())
      El::RuntimeError("Objectives not found in input files.");
    if(sdp.normalization.empty())
      El::RuntimeError("Normalization not found in input files.");

    for(const auto &matrix : sdp.matrices)
      {
        for(auto &pole : matrix.damped_rational.poles)
          {
            if(pole > 0)
              {
                throw std::runtime_error(
                  "All poles must be negative, but found '" + to_string(pole)
                  + "'");
              }
          }
      }
  }
}

PMWP_SDP::PMWP_SDP(const std::filesystem::path &input_path)
    : PMWP_SDP(std::vector{input_path})
{}

PMWP_SDP::PMWP_SDP(const std::vector<std::filesystem::path> &input_paths)
{
  num_matrices = 0;

  for(const auto &file : collect_files_expanding_nsv(input_paths))
    {
      // Simple round-robin for matrices across all files
      // TODO optimize https://github.com/davidsd/sdpb/issues/150
      auto should_parse_matrix
        = [this](size_t matrix_index) {
            return (this->num_matrices + matrix_index) % El::mpi::Size()
                   == El::mpi::Rank();
          };
      PMWP_File_Parse_Result result(file, should_parse_matrix);

      // TODO throw error if objectives or normalization are defined in several files?
      if(!result.objective.empty())
        objective = std::move(result.objective);
      if(!result.normalization.empty())
        normalization = std::move(result.normalization);

      for(auto &[index_in_file, matrix] : result.parsed_matrices)
        {
          size_t global_index = num_matrices + index_in_file;
          matrix_index_local_to_global.push_back(global_index);
          matrices.emplace_back(std::move(matrix));
        }

      num_matrices += result.num_matrices;
    }

  validate(*this);
}