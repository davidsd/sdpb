#include "PMP_File_Parse_Result.hxx"
#include "pmp_read.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdpb_util/block_mapping/compute_block_grid_mapping.hxx"
#include "sdpb_util/block_mapping/create_mpi_block_mapping_groups.hxx"

namespace fs = std::filesystem;

namespace
{
  struct Mapping
  {
    MPI_Comm_Wrapper mpi_comm;
    MPI_Group_Wrapper mpi_group;
    std::vector<size_t> file_indices;

    explicit Mapping(const std::vector<fs::path> &input_files)
    {
      const size_t num_files = input_files.size();

      // Distribute files among processes according to file size.
      std::vector<Block_Cost> file_costs;
      file_costs.reserve(num_files);
      for(size_t i = 0; i < num_files; ++i)
        {
          file_costs.emplace_back(fs::file_size(input_files.at(i)), i);
        }

      // We do not need shared memory windows, fast communication inside mpi_group etc.,
      // so we distribute files as if all processes were on a single node.
      const auto mapping
        = compute_block_grid_mapping(El::mpi::Size(), 1, file_costs);
      assert(mapping.size() == 1);

      create_mpi_block_mapping_groups(mapping, El::mpi::COMM_WORLD, 0,
                                      mpi_group.value, mpi_comm.value,
                                      file_indices);
    }
  };

  // Broadcast vector from the first rank for which vector is not empty.
  template <class T>
  void synchronize_vector(std::vector<T> &vec,
                          const El::mpi::Comm &comm = El::mpi::COMM_WORLD)
  {
    // Choose the first rank that has data
    int source_rank = vec.empty() ? comm.Size() : comm.Rank();
    source_rank = El::mpi::AllReduce(source_rank, El::mpi::MIN, comm);

    // vector is empty on all ranks
    if(source_rank == comm.Size())
      return;

    size_t size = vec.size();
    El::mpi::Broadcast(size, source_rank, comm);

    vec.resize(size);
    El::mpi::Broadcast(vec.data(), size, source_rank, comm);
  }
}

Polynomial_Matrix_Program
read_polynomial_matrix_program(const fs::path &input_file)
{
  return read_polynomial_matrix_program(std::vector{input_file});
}

// Read Polynomal Matrix Program in one of the supported formats.
//
// File reading is parallelized as follows:
// 1. Input files are distributed among (groups of) processes
//    via compute_block_grid_mapping().
// 2. Within each group, matrices from input files are distributed
//    among processes in a round-robin way.
//    Each process reads its matrices and skips the rest.
// TODO: since IO is often a bottleneck, we can reduce it:
//   root of each group copies file content to a shared memory window,
//   and other processes read it.
Polynomial_Matrix_Program
read_polynomial_matrix_program(const std::vector<fs::path> &input_files)
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

  const auto all_files = collect_files_expanding_nsv(input_files);
  const size_t num_files = all_files.size();

  // Parse files

  std::map<size_t, PMP_File_Parse_Result> parse_results;
  {
    // Determine which processes will parse which files
    const Mapping mapping(all_files);
    // Cumulative number of matrices in all files parsed
    // by the current group of processes
    size_t num_matrices_in_group = 0;
    for(const auto file_index : mapping.file_indices)
      {
        const auto &file = all_files.at(file_index);
        // Simple round-robin for matrices across files in a given group
        auto should_parse_matrix
          = [&mapping, &num_matrices_in_group](size_t matrix_index_in_file) {
              return (num_matrices_in_group + matrix_index_in_file)
                       % mapping.mpi_comm.value.Size()
                     == mapping.mpi_comm.value.Rank();
            };
        // TODO set also bool should_parse_objective and should_parse_normalization
        // to (mapping.mpi_comm.value.Rank() == 0)

        auto file_parse_result
          = PMP_File_Parse_Result::read(file, should_parse_matrix);

        num_matrices_in_group += file_parse_result.num_matrices;

        {
          auto [it, inserted]
            = parse_results.insert({file_index, std::move(file_parse_result)});
          if(!inserted)
            El::LogicError("Parsing result already insterted: ",
                           all_files.at(file_index));
        }
      }
  }

  // Total number of matrices in previous files.
  // Used to calculate matrix_index_local_to_global (= index_local + offset)
  std::vector<size_t> matrix_index_offset_per_file(num_files, 0);

  // Synhronize number of matrices in each file
  {
    std::vector<size_t> num_matrices_per_file(num_files, 0);
    for(auto &[file_index, parse_result] : parse_results)
      {
        num_matrices_per_file.at(file_index) = parse_result.num_matrices;
      }
    El::mpi::AllReduce(num_matrices_per_file.data(), num_files, El::mpi::MAX,
                       El::mpi::COMM_WORLD);

    matrix_index_offset_per_file.at(0) = 0;
    for(size_t file_index = 1; file_index < num_files; ++file_index)
      {
        matrix_index_offset_per_file.at(file_index)
          = matrix_index_offset_per_file.at(file_index - 1)
            + num_matrices_per_file.at(file_index - 1);
      }

    num_matrices
      = std::reduce(num_matrices_per_file.begin(), num_matrices_per_file.end(),
                    static_cast<size_t>(0));
  }

  // Get objective, normalization and polynomial vector matrices
  // + calculate block indices
  // NB: after this loop, data it moved from parse_results, don't use it!
  for(auto &&[file_index, parse_result] : parse_results)
    {
      if(!parse_result.objective.empty())
        {
          if(!objective.empty())
            El::RuntimeError("Duplicate objective at ",
                             all_files.at(file_index));
          objective = std::move(parse_result.objective);
        }
      if(!parse_result.normalization.empty())
        {
          if(!normalization.empty())
            El::RuntimeError("Duplicate normalization at ",
                             all_files.at(file_index));
          normalization = std::move(parse_result.normalization);
        }
      for(auto &[index_in_file, matrix] : parse_result.parsed_matrices)
        {
          size_t global_index
            = matrix_index_offset_per_file.at(file_index) + index_in_file;
          matrices.emplace_back(std::move(matrix));
          matrix_index_local_to_global.push_back(global_index);
        }
    }

  // TODO check also that objective/normalization exist only in a single file
  synchronize_vector(objective);
  synchronize_vector(normalization);

  if(objective.empty())
    El::RuntimeError("objective not found in input files");

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