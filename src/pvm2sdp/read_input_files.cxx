#include "sdp_convert/Dual_Constraint_Group.hxx"
#include "sdp_read/polynomial_vector_matrices/PVM_SDP.hxx"

namespace fs = std::filesystem;

void read_input_files(
  const std::vector<fs::path> &input_files, El::BigFloat &objective_const,
  std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  size_t &num_processed)
{
  PVM_SDP sdp(input_files);
  // b_0
  objective_const = sdp.objective.at(0);
  // b_1 ..b_N
  dual_objectives_b = std::vector<El::BigFloat>(sdp.objective.begin() + 1,
                                                sdp.objective.end());

  for(size_t i = 0; i < sdp.matrices.size(); ++i)
    {
      size_t block_index = sdp.matrix_index_local_to_global.at(i);
      auto &matrix = sdp.matrices.at(i);
      dual_constraint_groups.emplace_back(block_index, matrix);
      // Free memory
      matrix.clear();
    }
  num_processed = sdp.num_matrices;
}
