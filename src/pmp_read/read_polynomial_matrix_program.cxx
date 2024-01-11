#include "PMP_File_Parse_Result.hxx"
#include "pmp_read.hxx"

Polynomial_Matrix_Program
read_polynomial_matrix_program(const std::filesystem::path &input_file)
{
  return read_polynomial_matrix_program(std::vector{input_file});
}
// Read Polynomal Matrix Program
// in one of the supported formats
Polynomial_Matrix_Program read_polynomial_matrix_program(
  const std::vector<std::filesystem::path> &input_files)
{
  std::vector<El::BigFloat> objective;
  std::vector<El::BigFloat> normalization;
  // Total number of PVM matrices
  size_t num_matrices = 0;
  // In case of several processes,
  // each process owns only some matrices.
  std::vector<Polynomial_Vector_Matrix> matrices;
  // global index of matrices[i], lies in [0..num_matrices)
  std::vector<size_t> matrix_index_local_to_global;

  for(const auto &file : collect_files_expanding_nsv(input_files))
    {
      // Simple round-robin for matrices across all files
      // TODO optimize https://github.com/davidsd/sdpb/issues/150
      auto should_parse_matrix = [&num_matrices](size_t matrix_index) {
        return (num_matrices + matrix_index) % El::mpi::Size()
               == El::mpi::Rank();
      };
      PMP_File_Parse_Result file_parse_result(file, should_parse_matrix);

      // TODO throw error if objectives or normalization are defined in several files?
      if(!file_parse_result.objective.empty())
        objective = std::move(file_parse_result.objective);
      if(!file_parse_result.normalization.empty())
        normalization = std::move(file_parse_result.normalization);

      for(auto &[index_in_file, matrix] : file_parse_result.parsed_matrices)
        {
          size_t global_index = num_matrices + index_in_file;
          matrices.emplace_back(std::move(matrix));
          matrix_index_local_to_global.push_back(global_index);
        }

      num_matrices += file_parse_result.num_matrices;
    }

  if(normalization.empty())
    {
      // default normalization [1,0,0..0]
      normalization.resize(objective.size(), 0);
      normalization.at(0) = 1;
    }

  return Polynomial_Matrix_Program{
    std::move(objective), std::move(normalization), num_matrices,
    std::move(matrices), std::move(matrix_index_local_to_global)};
}