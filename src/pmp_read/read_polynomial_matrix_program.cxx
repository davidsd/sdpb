#include "PMP_File_Parse_Result.hxx"
#include "pmp_read.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdpb_util/assert.hxx"
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

    Mapping(const Environment &env, const std::vector<fs::path> &input_files)
    {
      const size_t num_files = input_files.size();

      // Distribute files among processes according to file size.
      std::vector<Block_Cost> file_costs;
      file_costs.reserve(num_files);
      for(size_t i = 0; i < num_files; ++i)
        {
          file_costs.emplace_back(fs::file_size(input_files.at(i)), i);
        }

      const auto node_comm = env.comm_shared_mem;
      ASSERT_EQUAL(node_comm.Size() * env.num_nodes(), El::mpi::Size(),
                   DEBUG_STRING(node_comm.Size()),
                   DEBUG_STRING(env.num_nodes()),
                   "Each node should have the same number of processes.");
      // Reading file from different nodes may be slower,
      // so we distribute files across nodes and processes
      // in such a way that each file is assigned to a single node.
      // NB: in case of a single input file, only first node will process it.
      // It can be slower than reading it from all ranks.
      // However, reading the same file from too many ranks looks non-optimal too.
      // In that case user should split input to several files.
      const auto mapping = compute_block_grid_mapping(
        node_comm.Size(), env.num_nodes(), file_costs);
      ASSERT_EQUAL(mapping.size(), env.num_nodes());

      create_mpi_block_mapping_groups(mapping, node_comm, env.node_index(),
                                      mpi_group.value, mpi_comm.value,
                                      file_indices);
    }
  };

  // Check that vec.value() is the same for all ranks where vec.has_value()
  // and broadcast it to other ranks.
  template <class T>
  [[nodiscard]] std::optional<std::vector<T>>
  check_and_broadcast_vector(const std::optional<std::vector<T>> &vec,
                             // Parameters for debug output:
                             const std::string &vector_name,
                             const std::optional<size_t> &file_index,
                             const std::vector<std::filesystem::path> &files,
                             const El::mpi::Comm &comm = El::mpi::COMM_WORLD)
  {
    if(comm.Size() == 1)
      return vec;

    // Choose the first rank that has data
    int source_rank = vec.has_value() ? comm.Rank() : comm.Size();
    source_rank = El::mpi::AllReduce(source_rank, El::mpi::MIN, comm);

    // vector is missing on all ranks
    if(source_rank == comm.Size())
      return vec;

    std::vector<T> result = vec.value_or(std::vector<T>{});
    // broadcast vector
    {
      size_t size = result.size();
      El::mpi::Broadcast(size, source_rank, comm);

      result.resize(size);
      El::mpi::Broadcast(result.data(), size, source_rank, comm);
    }

    // Check that all non-empty vectors were the same
    {
      size_t source_file_index = file_index.value_or(files.size());
      El::mpi::Broadcast(source_file_index, source_rank, comm);

      if(vec.has_value())
        {
          ASSERT(result == vec.value(), "Found different ", vector_name,
                 " vectors in input files:\n\t", files.at(file_index.value()),
                 "\n\t", files.at(source_file_index));
        }
    }
    return result;
  }
}

Polynomial_Matrix_Program
read_polynomial_matrix_program(const Environment &env,
                               const fs::path &input_file,
                               const Verbosity &verbosity, Timers &timers)
{
  return read_polynomial_matrix_program(env, std::vector{input_file},
                                        verbosity, timers);
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
read_polynomial_matrix_program(const Environment &env,
                               const std::vector<fs::path> &input_files,
                               const Verbosity &verbosity, Timers &timers)
{
  Scoped_Timer timer(timers, "read_pmp");

  std::optional<std::vector<El::BigFloat>> objective;
  std::optional<std::vector<El::BigFloat>> normalization;
  // Total number of PVM matrices
  size_t num_matrices = 0;
  // In case of several processes,
  // each process owns only some matrices.
  std::vector<Polynomial_Vector_Matrix> matrices;
  // global index of matrices[i], lies in [0..num_matrices)
  std::vector<size_t> matrix_index_local_to_global;
  // input path for each of the matrices
  std::vector<std::filesystem::path> block_paths;

  const auto all_files = collect_files_expanding_nsv(input_files);
  const size_t num_files = all_files.size();

  // Parse files

  std::map<size_t, PMP_File_Parse_Result> parse_results;
  {
    Scoped_Timer parse_timer(timers, "parse");
    // Determine which processes will parse which files
    const Mapping mapping(env, all_files);
    // Cumulative number of matrices in all files parsed
    // by the current group of processes
    size_t num_matrices_in_group = 0;
    for(const auto file_index : mapping.file_indices)
      {
        const auto &file = all_files.at(file_index);
        Scoped_Timer parse_file_timer(timers,
                                      "file_" + std::to_string(file_index)
                                        + "=" + file.filename().string());
        if(verbosity >= Verbosity::trace)
          {
            El::Output("rank=", El::mpi::Rank(), " read ", file);
          }

        // Simple round-robin for matrices across files in a given group
        auto should_parse_matrix
          = [&mapping, &num_matrices_in_group](size_t matrix_index_in_file) {
              return (num_matrices_in_group + matrix_index_in_file)
                       % mapping.mpi_comm.value.Size()
                     == mapping.mpi_comm.value.Rank();
            };

        bool should_parse_objective = mapping.mpi_comm.value.Rank() == 0;
        bool should_parse_normalization = mapping.mpi_comm.value.Rank() == 0;

        auto file_parse_result = PMP_File_Parse_Result::read(
          file, should_parse_objective, should_parse_normalization,
          should_parse_matrix);

        num_matrices_in_group += file_parse_result.num_matrices;

        {
          auto [it, inserted]
            = parse_results.insert({file_index, std::move(file_parse_result)});
          ASSERT(inserted, "Duplicate parsing result for ",
                 all_files.at(file_index));
        }
      }
  }

  // Total number of matrices in previous files.
  // Used to calculate matrix_index_local_to_global (= index_local + offset)
  std::vector<size_t> matrix_index_offset_per_file(num_files, 0);

  // Synhronize number of matrices in each file
  {
    Scoped_Timer sync_num_matrices_timer(timers, "sync_num_matrices");

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
      = std::accumulate(num_matrices_per_file.begin(),
                        num_matrices_per_file.end(), static_cast<size_t>(0));
  }

  // Get objective, normalization and polynomial vector matrices
  // + calculate block indices
  // NB: after this loop, data it moved from parse_results, don't use it!
  std::optional<size_t> objective_file_index;
  std::optional<size_t> normalization_file_index;
  for(auto &&[file_index, parse_result] : parse_results)
    {
      if(parse_result.objective.has_value())
        {
          if(!objective.has_value())
            {
              objective = std::move(parse_result.objective);
              objective_file_index = file_index;
            }
          else
            {
              ASSERT(objective == parse_result.objective,
                     "Found different objective vectors in input files:\n\t",
                     all_files.at(objective_file_index.value()), "\n\t",
                     all_files.at(file_index));
            }
        }
      if(parse_result.normalization.has_value())
        {
          if(!normalization.has_value())
            {
              normalization = std::move(parse_result.normalization);
              normalization_file_index = file_index;
            }
          else
            {
              ASSERT(
                normalization == parse_result.normalization,
                "Found different normalization vectors in input files:\n\t",
                all_files.at(normalization_file_index.value()), "\n\t",
                all_files.at(file_index));
            }
        }
      for(auto &[index_in_file, matrix] : parse_result.parsed_matrices)
        {
          size_t global_index
            = matrix_index_offset_per_file.at(file_index) + index_in_file;
          matrices.emplace_back(std::move(matrix));
          matrix_index_local_to_global.push_back(global_index);
          block_paths.push_back(all_files.at(file_index));
        }
    }

  {
    Scoped_Timer sync_objective_timer(timers, "sync_objective_normalization");
    // TODO we can store objective on one rank, no need to synchronize it
    objective = check_and_broadcast_vector(objective, "objective",
                                           objective_file_index, all_files);
    normalization = check_and_broadcast_vector(
      normalization, "normalization", normalization_file_index, all_files);
  }

  ASSERT(objective.has_value(), "objective not found in input files");

  return Polynomial_Matrix_Program(
    std::move(objective.value()), std::move(normalization), num_matrices,
    std::move(matrices), std::move(matrix_index_local_to_global),
    std::move(block_paths));
}