#include "Approx_Parameters.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/initialize_bigint_syrk_context.hxx"
#include "sdpb_util/memory_estimates.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <filesystem>

namespace fs = std::filesystem;

// TODO: Have this be part of sdp_solve.hxx
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L,
                            const Block_Info &block_info,
                            const std::string &name);
void compute_A_X_inv(
  const Block_Info &block_info, const Block_Diagonal_Matrix &X_cholesky,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  std::array<std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>,
             2> &A_X_inv);

void compute_A_Y(
  const Block_Info &block_info, const Block_Diagonal_Matrix &Y,
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
  // Estimate how many BigFloats will be allocated by SDPB on the current node,
  // (including what's already allocated, e.g. SDP)
  size_t get_required_nonshared_memory_per_node_bytes(
    const Environment &env, const Block_Info &block_info, const SDP &sdp,
    const Block_Diagonal_Matrix &X, const Verbosity verbosity)
  {
    const auto &node_comm = env.comm_shared_mem;

    // X, Y, X_cholesky, Y_cholesky, primal_residues
    const size_t X_size
      = El::mpi::Reduce(get_matrix_size_local(X), 0, node_comm);

    // Bilinear pairing blocks - A_X_inv, A_Y
    const size_t A_X_inv_size
      = El::mpi::Reduce(get_A_X_size_local(block_info, sdp), 0, node_comm);

    // schur_complement, schur_complement_cholesky
    const size_t schur_complement_size = El::mpi::Reduce(
      get_schur_complement_size_local(block_info), 0, node_comm);

    // #(B) = PxN
    // free_var_matrix, schur_off_diagonal
    const size_t B_size = El::mpi::Reduce(get_B_size_local(sdp), 0, node_comm);

    // #Q = NxN, distributed over all nodes.
    const size_t Q_size = El::mpi::Reduce(get_Q_size_local(sdp), 0, node_comm);

    // Memory for new SDP created in quadratic_approximate_objectives()
    const size_t SDP_size
      = El::mpi::Reduce(get_SDP_size_local(sdp), 0, node_comm);

    // We will use only result on rank=0
    if(node_comm.Rank() != 0)
      return 0;

    // Calculate mem_required_size
    size_t mem_required_size = 0;

    // SDP struct
    mem_required_size += SDP_size;

    // X, Y, X_cholesky, Y_cholesky, primal_residues
    mem_required_size += 5 * X_size;

    // A_X_inv and A_Y
    mem_required_size += 2 * A_X_inv_size;

    // schur_complement_cholesky
    mem_required_size += schur_complement_size;

    // Add either schur_complement from initialize_schur_complement_solver(),
    // or SDP from quadratic_approximate_objectives()
    // (they do not coexist, thus we choose maximum size instead of adding both)
    mem_required_size += std::max(schur_complement_size, SDP_size);

    // schur_off_diagonal = L^{-1} B takes the same size as B
    // Allocated in initialize_schur_off_diagonal()
    // (B = sdp.free_var_matrix is already allocated)
    mem_required_size += B_size;
    // Q = NxN
    mem_required_size += Q_size;

    // initial_node_mem_used() is RAM allocated at SDPB start.
    // This could be important: e.g. on 128 cores (Expanse HPC) it is ~26GB
    const size_t mem_required_bytes
      = env.initial_node_mem_used() + mem_required_size * bigfloat_bytes();

    if(verbosity >= Verbosity::debug)
      {
        std::ostringstream ss;
        El::BuildStream(
          ss, "node=", env.node_index(),
          " matrix sizes and memory estimates: ", "\n\t#(SDP) = ", SDP_size,
          "\n\t#(X) = ", X_size, "\n\t#(A_X_inv) = ", A_X_inv_size,
          "\n\t#(schur_complement) = ", schur_complement_size,
          "\n\t#(B) = ", B_size, "\n\t#(Q) = ", Q_size,
          "\n\tBigFloat size: ", pretty_print_bytes(bigfloat_bytes()),
          "\n\tTotal BigFloats to be allocated: ", mem_required_size,
          " elements = ",
          pretty_print_bytes(mem_required_size * bigfloat_bytes()),
          "\n\tInitial MemUsed (at SDPB start) = ",
          pretty_print_bytes(env.initial_node_mem_used()),
          "\n\tTotal non-shared memory estimate: ",
          pretty_print_bytes(mem_required_bytes, true));
        El::Output(ss.str());
      }

    return mem_required_bytes;
  }

  size_t
  get_max_shared_memory_bytes(const size_t default_max_shared_memory_bytes,
                              const Environment &env,
                              const Block_Info &block_info, const SDP &sdp,
                              const Block_Diagonal_Matrix &X,
                              const Verbosity verbosity)
  {
    // If user sets --maxSharedMemory limit manually, we use it.
    // Otherwise, we calculate the limit automatically.
    if(default_max_shared_memory_bytes != 0)
      return default_max_shared_memory_bytes;
    const size_t nonshared_memory_required_per_node_bytes
      = get_required_nonshared_memory_per_node_bytes(env, block_info, sdp, X,
                                                     verbosity);
    return get_max_shared_memory_bytes(
      nonshared_memory_required_per_node_bytes, env, verbosity);
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
      Block_Diagonal_Matrix X(block_info.psd_matrix_block_sizes(),
                              block_info.block_indices,
                              block_info.num_points.size(), grid),
        Y(X);
      for(size_t block = 0; block != block_info.block_indices.size(); ++block)
        {
          size_t block_index(block_info.block_indices.at(block));
          for(size_t psd_block(0); psd_block < 2; ++psd_block)
            {
              // Constant constraints have empty odd parity blocks, so we do
              // not need to load them.
              if(X.blocks.at(2 * block + psd_block).Height() != 0)
                {
                  const size_t psd_index(2 * block_index + psd_block);
                  read_text_block(X.blocks.at(2 * block + psd_block),
                                  solution_dir, "X_matrix_", psd_index);
                  read_text_block(Y.blocks.at(2 * block + psd_block),
                                  solution_dir, "Y_matrix_", psd_index);
                }
            }
        }

      std::array<
        std::vector<std::vector<std::vector<El::DistMatrix<El::BigFloat>>>>, 2>
        A_X_inv, A_Y;

      Block_Diagonal_Matrix X_cholesky(X);
      cholesky_decomposition(X, X_cholesky, block_info, "X");
      compute_A_X_inv(block_info, X_cholesky, sdp.bases_blocks, A_X_inv);
      compute_A_Y(block_info, Y, sdp.bases_blocks, A_Y);

      Timers timers;
      El::Matrix<int32_t> block_timings_ms(block_info.dimensions.size(), 1);
      El::Zero(block_timings_ms);

      const auto verbosity = parameters.verbosity;
      const auto max_shared_memory_bytes
        = get_max_shared_memory_bytes(parameters.max_shared_memory_bytes, env,
                                      block_info, sdp, X, verbosity);
      auto bigint_syrk_context = initialize_bigint_syrk_context(
        env, block_info, sdp, max_shared_memory_bytes, verbosity);

      initialize_schur_complement_solver(
        env, block_info, sdp, A_X_inv, A_Y, grid, schur_complement_cholesky,
        schur_off_diagonal, bigint_syrk_context, Q, timers, block_timings_ms,
        verbosity);
    }
}
