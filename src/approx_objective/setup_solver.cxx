#include "Approx_Parameters.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "sdp_solve/SDP_Solver/run/get_syrk_Q_config.hxx"
#include "sdp_solve/memory_estimates.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <filesystem>

namespace fs = std::filesystem;

// TODO: Have this be part of sdp_solve.hxx
void cholesky_decomposition(const Paired_Block_Diagonal_Matrix &A,
                            Paired_Block_Diagonal_Matrix &L,
                            const Block_Info &block_info,
                            const std::string &name);
void compute_A_X_inv(
  const Block_Info &block_info, const Paired_Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_X_inv);

void compute_A_Y(
  const Block_Info &block_info, const Paired_Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_Y);

void initialize_schur_complement_solver(
  const Environment &env, const Block_Info &block_info, const SDP &sdp,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_X_inv,
  const std::array<
    std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
    &A_Y,
  const El::Grid &block_grid, Block_Diagonal_Matrix &schur_complement_cholesky,
  Block_Matrix &schur_off_diagonal,
  BigInt_Shared_Memory_Syrk_Context &bigint_syrk_context,
  El::DistMatrix<El::BigFloat> &Q, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms, Verbosity verbosity);

namespace
{
  BigInt_Shared_Memory_Syrk_Context create_syrk_Q_context_and_print_memory(
    const Environment &env, const Block_Info &block_info, const SDP &sdp,
    const Paired_Block_Diagonal_Matrix &X, const Approx_Parameters &parameters)
  {
    // Estimate how many BigFloats will be allocated by SDPB on the current node,
    // (including what's already allocated, e.g. SDP)
    const auto &node_comm = env.comm_shared_mem;
    const auto node_reduce = [&node_comm](const size_t value) -> size_t {
      return El::mpi::AllReduce(value, node_comm);
    };

    // X, Y, X_cholesky, Y_cholesky, primal_residues
    const size_t X_size = node_reduce(get_matrix_size_local(X));
    const size_t X_bytes = node_reduce(get_allocated_bytes(X));

    // Bilinear pairing blocks - A_X_inv, A_Y
    const size_t A_X_inv_size
      = El::mpi::Reduce(get_A_X_size_local(block_info, sdp), 0, node_comm);
    const size_t A_X_inv_bytes = A_X_inv_size * bigfloat_bytes();

    // schur_complement, schur_complement_cholesky
    const size_t schur_complement_size
      = node_reduce(get_schur_complement_size_local(block_info));
    const size_t schur_complement_bytes
      = schur_complement_size * bigfloat_bytes();

    // #(B) = PxN
    // free_var_matrix, schur_off_diagonal
    const size_t B_size = node_reduce(get_B_size_local(sdp));
    const size_t B_bytes
      = node_reduce(get_allocated_bytes(sdp.free_var_matrix));

    // #Q = NxN, distributed over all nodes.
    const size_t Q_size = node_reduce(get_Q_size_local(sdp));
    const size_t Q_bytes = Q_size * bigfloat_bytes();

    // Memory for new SDP created in quadratic_approximate_objectives()
    const size_t SDP_bytes = node_reduce(get_allocated_bytes(sdp));

    size_t solver_bytes = 0;

    const auto mem_required_bytes = [&] { return SDP_bytes + solver_bytes; };

    // X, Y, X_cholesky, Y_cholesky, primal_residues
    solver_bytes += 5 * X_bytes;

    // A_X_inv and A_Y
    solver_bytes += 2 * A_X_inv_bytes;

    // schur_complement_cholesky
    solver_bytes += schur_complement_bytes;

    // schur_off_diagonal = L^{-1} B takes the same size as B
    // Allocated in initialize_schur_off_diagonal()
    // (B = sdp.free_var_matrix is already allocated)
    solver_bytes += B_bytes;
    // Q = NxN
    solver_bytes += Q_bytes;

    const auto cfg = get_syrk_Q_config(
      block_info, sdp, parameters.max_memory, parameters.max_shared_memory,
      mem_required_bytes() + std::max(schur_complement_bytes, SDP_bytes));

    // Add either schur_complement from initialize_schur_complement_solver(),
    // or SDP from quadratic_approximate_objectives()
    // (they do not coexist, thus we choose maximum size instead of adding both)
    solver_bytes
      += std::max(schur_complement_bytes + cfg.node_local_bytes(), SDP_bytes);
    // Shared memory windows are allocated in the beginning, so they should be always included.
    solver_bytes += cfg.node_shmem_bytes();

    if(node_comm.Rank() == 0 && parameters.verbosity >= Verbosity::debug)
      {
        // Print memory estimates

        std::vector<std::pair<std::string, size_t>> num_elements_per_category{
          {"X", X_size},
          {"A_X_inv", A_X_inv_size},
          {"schur_complement", schur_complement_size},
          {"B", B_size},
          {"Q", Q_size},
        };

        std::vector<std::pair<std::string, size_t>> bytes_per_category{
          {"Initial MemAvailable (at SDPB start)",
           env.initial_node_mem_available()},
          {"BigFloat size", bigfloat_bytes()},
          {"Total SDPB memory estimate", mem_required_bytes()},
          {"Shared memory estimate", cfg.node_shmem_bytes()},
          {"\tSolver", solver_bytes},
          {"\t\tSDP", SDP_bytes},
          {"\t\tinitialize_schur_complement_solver()",
           schur_complement_bytes + cfg.node_total_bytes()},
          {"\t\t\tschur_complement", schur_complement_bytes},
          {"\t\t\tQ", Q_bytes},
          {"\t\t\t\tshared memory", cfg.node_shmem_bytes()},
        };

        std::ostringstream ss;
        El::BuildStream(ss, "node=", env.node_index(),
                        " matrix sizes and memory estimates: ");

        for(const auto &[name, size] : num_elements_per_category)
          {
            El::BuildStream(ss, "\n\t#(", name, ") = ", size, " elements");
          }
        for(const auto &[name, bytes] : bytes_per_category)
          {
            El::BuildStream(ss, "\n\t", name, ": ",
                            pretty_print_bytes(bytes, true));
          }
        El::BuildStream(
          ss, "\n\tShared memory configuration: ",
          "\n\t\tsyrk_Q() input split factor: ", cfg.input_split_factor,
          "\n\t\tsyrk_Q() output split factor: ", cfg.output_split_factor);
        El::Output(ss.str());
      }

    return BigInt_Shared_Memory_Syrk_Context(
      cfg, block_info.node_group_index(), block_info.block_indices,
      parameters.verbosity);
  }
}

void setup_solver(const Environment &env, const Block_Info &block_info,
                  const El::Grid &grid, const SDP &sdp,
                  const Approx_Parameters &parameters,
                  Block_Diagonal_Matrix &schur_complement_cholesky,
                  Block_Matrix &schur_off_diagonal,
                  El::DistMatrix<El::BigFloat> &Q)
{
  const auto &solution_dir = parameters.solution_dir;
  if(fs::exists(solution_dir / "Q_cholesky.txt"))
    {
      for(size_t block = 0; block != block_info.block_indices.size(); ++block)
        {
          size_t block_index(block_info.block_indices.at(block));
          read_text_block(schur_complement_cholesky.blocks.at(block),
                          solution_dir, "schur_complement_cholesky_",
                          block_index);
          read_text_block(schur_off_diagonal.blocks.at(block), solution_dir,
                          "schur_off_diagonal_", block_index);
        }
      read_text_block(Q, solution_dir / "Q_cholesky.txt");
    }
  else
    {
      Paired_Block_Diagonal_Matrix X(block_info.psd_matrix_block_sizes(),
                                     block_info.block_indices, grid),
        Y(X);
      for(size_t block = 0; block != block_info.block_indices.size(); ++block)
        {
          size_t block_index(block_info.block_indices.at(block));
          for(size_t psd_block(0); psd_block < 2; ++psd_block)
            {
              const size_t psd_index = 2 * block_index + psd_block;
              auto &X_block = X.get_block(block, psd_block);
              auto &Y_block = Y.get_block(block, psd_block);
              // Constant constraints have empty odd parity blocks, so we do
              // not need to load them.
              if(X_block.Height() != 0)
                {
                  read_text_block(X_block, solution_dir, "X_matrix_",
                                  psd_index);
                  read_text_block(Y_block, solution_dir, "Y_matrix_",
                                  psd_index);
                }
            }
        }

      std::array<
        std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
        A_X_inv, A_Y;

      auto X_cholesky = X;
      cholesky_decomposition(X, X_cholesky, block_info.block_indices, "X");
      compute_A_X_inv(block_info, X_cholesky, sdp.bases_blocks, A_X_inv);
      compute_A_Y(block_info, Y, sdp.bases_blocks, A_Y);

      Timers timers;
      El::Matrix<int32_t> block_timings_ms(block_info.dimensions.size(), 1);
      El::Zero(block_timings_ms);

      const auto verbosity = parameters.verbosity;
      const auto max_memory = parameters.max_shared_memory;
      auto bigint_syrk_context = create_syrk_Q_context_and_print_memory(
        env, block_info, sdp, X, parameters);

      initialize_schur_complement_solver(
        env, block_info, sdp, A_X_inv, A_Y, grid, schur_complement_cholesky,
        schur_off_diagonal, bigint_syrk_context, Q, timers, block_timings_ms,
        verbosity);
    }
}
