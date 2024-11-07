#include "../BigInt_Shared_Memory_Syrk_Context.hxx"
#include "../fmpz/fmpz_mul_blas_util.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <El.hpp>

#include <cblas.h>

// code adopted from flint mul_blas.c
namespace
{
  [[nodiscard]] El::Int sum(const std::vector<El::Int> &block_heights)
  {
    return std::accumulate(block_heights.begin(), block_heights.end(), 0);
  }

  template <class T = double>
  [[nodiscard]] size_t
  window_size_bytes(const size_t height, const size_t width,
                    const size_t num_primes)
  {
    return height * width * num_primes * sizeof(T);
  }

  template <class T>
  [[nodiscard]] size_t
  window_size_bytes(const Residue_Matrices_Window<T> &window)
  {
    return window_size_bytes<T>(window.height, window.width,
                                window.num_primes);
  }

  // Divide and round to the nearest integer above
  [[nodiscard]] size_t div_ceil(const size_t x, const size_t y)
  {
    return std::ceil(static_cast<double>(x) / static_cast<double>(y));
  }

  // Estimate total MPI buffer sizes on a node for restore_and_reduce()
  size_t
  get_reduce_scatter_buffer_bytes(const El::mpi::Comm &shared_memory_comm,
                                  const size_t window_width,
                                  const size_t split_factor)
  {
    ASSERT(El::mpi::Size() % shared_memory_comm.Size() == 0,
           "All nodes should have the same number of ranks! ",
           DEBUG_STRING(El::mpi::Size()),
           DEBUG_STRING(shared_memory_comm.Size()));
    auto num_nodes = El::mpi::Size() / shared_memory_comm.Size();

    // no reduce-scatter => no RAM needed
    if(num_nodes == 1)
      return 0;

    // Total number of elements in Q submatrix
    // If split factor == 1, we need only half of the matrix.
    // Otherwise, we need the whole matrix for off-diagonal blocks
    const size_t num_elements = split_factor == 1
                                  ? window_width * (window_width + 1) / 2
                                  : window_width * window_width;
    const size_t bigfloat_bytes = El::BigFloat(1.1).SerializedSize();

    // Each rank needs ~(num_elements/num_ranks) buffer for MPI_Send
    // and the same for MPI_Recv
    return 2 * div_ceil(num_elements, num_nodes) * bigfloat_bytes;
  }

  [[nodiscard]] bool calculate_input_window_split(
    const std::vector<El::Int> &blocks_height_per_group,
    const size_t max_input_window_height_per_prime,
    std::vector<El::Int> &input_window_height_per_group_per_prime,
    size_t &input_window_split_factor)
  {
    const size_t num_groups = blocks_height_per_group.size();
    const size_t total_block_height_per_node = sum(blocks_height_per_group);

    if(max_input_window_height_per_prime < num_groups)
      return false;

    const size_t min_split_factor = div_ceil(
      total_block_height_per_node, max_input_window_height_per_prime);

    // Largest possible split factor is max(group heights)
    const size_t max_split_factor = std::accumulate(
      blocks_height_per_group.begin(), blocks_height_per_group.end(), 0,
      [](auto a, auto b) { return std::max(a, b); });

    for(input_window_split_factor = min_split_factor;
        input_window_split_factor <= max_split_factor;
        ++input_window_split_factor)
      {
        size_t window_height = 0;
        for(size_t i = 0; i < num_groups; ++i)
          {
            const auto curr_height = div_ceil(blocks_height_per_group.at(i),
                                              input_window_split_factor);
            input_window_height_per_group_per_prime.at(i) = curr_height;
            window_height += curr_height;
          }
        if(window_height <= max_input_window_height_per_prime)
          return true;
      }
    LOGIC_ERROR(
      "Failed to calculate input window splitting for split_factors from ",
      min_split_factor, " to ", max_split_factor, ". ",
      DEBUG_STRING(total_block_height_per_node),
      DEBUG_STRING(max_input_window_height_per_prime));
  }
}

BigInt_Shared_Memory_Syrk_Context::BigInt_Shared_Memory_Syrk_Context(
  const El::mpi::Comm &shared_memory_comm, const size_t group_index,
  const std::vector<int> &group_comm_sizes, const mp_bitcnt_t precision,
  size_t max_shared_memory_bytes,
  const std::vector<El::Int> &blocks_height_per_group, const int block_width,
  const std::vector<size_t> &block_index_local_to_global,
  const Verbosity verbosity,
  const std::function<Blas_Job_Schedule(
    Blas_Job::Kind kind, El::UpperOrLower uplo, size_t num_ranks,
    size_t num_primes, int output_height, int output_width,
    Verbosity _verbosity)> &create_job_schedule)
    : shared_memory_comm(shared_memory_comm),
      group_index(group_index),
      group_comm_sizes(group_comm_sizes),
      num_groups(group_comm_sizes.size()),
      total_block_height_per_node(sum(blocks_height_per_group)),
      comb(precision, precision, 1, total_block_height_per_node),
      verbosity(verbosity),
      block_index_local_to_global(block_index_local_to_global),
      create_blas_job_schedule_func(create_job_schedule)
{
  ASSERT_EQUAL(blocks_height_per_group.size(), num_groups);

  std::vector<int> input_window_height_per_group_per_prime(num_groups);
  size_t window_width;

  if(max_shared_memory_bytes == 0)
    max_shared_memory_bytes = std::numeric_limits<size_t>::max();

  size_t reduce_scatter_buffer_bytes;

  // Each extra split for output window leads to more reduce-scatter calls.
  // Thus, we try to find minimal output_window_split_factor
  // that allows to fit all shared windows into memory
  // (with potentially large input_window_split_factor, which is less critical for performance)
  // See details in ../Readme.md
  for(output_window_split_factor = 1;
      output_window_split_factor <= block_width; ++output_window_split_factor)
    {
      window_width = block_width / output_window_split_factor
                     + (block_width % output_window_split_factor == 0 ? 0 : 1);

      const auto output_window_bytes
        = window_size_bytes(window_width, window_width, comb.num_primes);

      reduce_scatter_buffer_bytes = get_reduce_scatter_buffer_bytes(
        shared_memory_comm, window_width, output_window_split_factor);

      // Try to find minimal input_window_split_factor
      size_t max_input_window_bytes;
      if(output_window_bytes + reduce_scatter_buffer_bytes
         >= max_shared_memory_bytes)
        {
          max_input_window_bytes = 0;
        }
      else
        {
          max_input_window_bytes = max_shared_memory_bytes
                                   - output_window_bytes
                                   - reduce_scatter_buffer_bytes;
        }

      // If output window is split, we need two (same-size) input windows
      // to calculate off-diagonal Q blocks
      if(output_window_split_factor > 1)
        max_input_window_bytes /= 2;

      const size_t residues_bytes_per_single_block_row
        = window_width * comb.num_primes * sizeof(double);
      const size_t max_input_window_height
        = max_input_window_bytes / residues_bytes_per_single_block_row;

      El::byte is_enough_memory = calculate_input_window_split(
        blocks_height_per_group, max_input_window_height,
        input_window_height_per_group_per_prime, input_window_split_factor);

      // output_window_split_factor should be the same for all nodes,
      // because we need global reduce-scatter for each output submatrix.
      // Thus, we require that is_enough_memory=true on all nodes.
      is_enough_memory = El::mpi::AllReduce(
        is_enough_memory, El::mpi::LOGICAL_AND, El::mpi::COMM_WORLD);
      if(is_enough_memory)
        break;
      // We tried maximal split factor, but still failed to fit:
      if(output_window_split_factor == block_width)
        {
          RUNTIME_ERROR(
            "Cannot allocate shared memory window for input block residues."
            " Required: at least ",
            pretty_print_bytes(residues_bytes_per_single_block_row, true),
            " per each MPI group, available: ",
            pretty_print_bytes(max_input_window_bytes), " in total. ",
            "\n\tmax_shared_memory = ",
            pretty_print_bytes(max_shared_memory_bytes, true),
            "\n\toutput_window = ",
            pretty_print_bytes(output_window_bytes, true),
            "\n\treduce_scatter_buffer = ",
            pretty_print_bytes(reduce_scatter_buffer_bytes, true), "\n\t",
            DEBUG_STRING(num_groups));
        }
    }

  ASSERT(output_window_split_factor > 0);
  ASSERT(input_window_split_factor > 0);

  // Print warnings for large split factors.
  if(shared_memory_comm.Rank() == 0)
    {
      const bool print_output_warning = output_window_split_factor > 1;
      const bool print_input_warning = input_window_split_factor > 10;
      if(print_output_warning || print_input_warning)
        {
          std::ostringstream ss;
          El::BuildStream(ss, "rank=", El::mpi::Rank(),
                          ": BigInt_Shared_Memory_Syrk_Context:\n");
          if(print_output_warning)
            {
              El::BuildStream(ss, "\tOutput window is split by a factor of ",
                              output_window_split_factor,
                              ", which may affect performance.\n");
            }
          if(print_input_warning)
            {
              El::BuildStream(ss,
                              "\tInput window is split by a large factor of ",
                              input_window_split_factor,
                              ", which may affect performance.\n");
            }
          El::BuildStream(
            ss,
            "\tConsider increasing available shared memory per node "
            "(--maxSharedMemory option).",
            "\n\tShared memory limit: ",
            pretty_print_bytes(max_shared_memory_bytes, true),
            "\n\tOutput window:", "\n\t\tactual size after splitting: ",
            pretty_print_bytes(
              window_size_bytes(window_width, window_width, comb.num_primes),
              true),
            "\n\t\toptimal size without splitting: ",
            pretty_print_bytes(
              window_size_bytes(block_width, block_width, comb.num_primes),
              true),
            "\n\tInput window:"
            "\n\t\tactual size after splitting: ",
            pretty_print_bytes(
              window_size_bytes(sum(input_window_height_per_group_per_prime),
                                window_width, comb.num_primes),
              true),
            "\n\t\toptimal size without splitting: ",
            pretty_print_bytes(window_size_bytes(total_block_height_per_node,
                                                 block_width, comb.num_primes),
                               true));
          PRINT_WARNING(ss.str());
        }
    }

  // Print sizes
  if(verbosity >= Verbosity::debug && shared_memory_comm.Rank() == 0)
    {
      std::ostringstream os;
      El::BuildStream(os, "create BigInt_Shared_Memory_Syrk_Context, rank=",
                      El::mpi::Rank(), "\n");
      El::BuildStream(os, "  Shared memory limit: ",
                      pretty_print_bytes(max_shared_memory_bytes, true), "\n");
      El::BuildStream(os, "  Number of primes: ", comb.num_primes, "\n");

      El::BuildStream(os, "  Blocks on the node:\n");
      El::BuildStream(os, "    Total height: ",
                      static_cast<const size_t>(total_block_height_per_node),
                      "\n");
      El::BuildStream(os, "    Width: ", block_width, "\n");
      El::BuildStream(
        os, "    Elements: ", total_block_height_per_node * block_width, "\n");
      // TODO for some reason, El::BuildStream() fails for std::vector
      // although we've included sdpb_util/ostream/ostream_vector.hxx
      os << "    Heights for each MPI group: " << blocks_height_per_group
         << "\n";

      El::BuildStream(os, "  Output residues window (Q):\n");
      El::BuildStream(
        os, "    Window size: ",
        pretty_print_bytes(
          window_size_bytes(window_width, window_width, comb.num_primes),
          true),
        "\n");
      El::BuildStream(os, "    Split factor: ", output_window_split_factor,
                      "\n");
      El::BuildStream(os, "    Height=Width (per prime): ", window_width,
                      "\n");
      El::BuildStream(os, "    Total elements (per prime): ",
                      window_width * window_width, "\n");

      auto input_window_height = sum(input_window_height_per_group_per_prime);

      El::BuildStream(os, "  Input residues window (P):\n");
      El::BuildStream(os, "    Number of windows: ",
                      output_window_split_factor == 1 ? 1 : 2, "\n");
      El::BuildStream(
        os, "    Window size: ",
        pretty_print_bytes(window_size_bytes(input_window_height, window_width,
                                             comb.num_primes),
                           true),
        "\n");
      El::BuildStream(os, "    Split factor: ", input_window_split_factor,
                      "\n");
      El::BuildStream(os, "    Height (per prime): ", input_window_height,
                      "\n");
      El::BuildStream(os, "    Width: ", window_width, "\n");
      os << "    Heights for each MPI group: "
         << input_window_height_per_group_per_prime << "\n";
      os << "    MPI group sizes:" << group_comm_sizes << "\n";

      if(reduce_scatter_buffer_bytes > 0)
        {
          El::BuildStream(
            os, "  Buffer size for reduce-scatter Q: ",
            pretty_print_bytes(reduce_scatter_buffer_bytes, true), "\n");
        }

      El::Output(os.str());
    }

  output_residues_window = std::make_unique<Residue_Matrices_Window<double>>(
    shared_memory_comm, comb.num_primes, window_width, window_width);

  input_grouped_block_residues_window_A
    = std::make_unique<Block_Residue_Matrices_Window<double>>(
      shared_memory_comm, comb.num_primes, num_groups,
      input_window_height_per_group_per_prime, window_width);

  // We need a second input window only to calculate off-diagonal blocks
  // of the output window.
  if(output_window_split_factor > 1)
    {
      input_grouped_block_residues_window_B
        = std::make_unique<Block_Residue_Matrices_Window<double>>(
          shared_memory_comm, comb.num_primes, num_groups,
          input_window_height_per_group_per_prime, window_width);
    }

  {
    // Check sizes
    auto total_bytes
      = window_size_bytes(*output_residues_window)
        + window_size_bytes(*input_grouped_block_residues_window_A);
    if(input_grouped_block_residues_window_B != nullptr)
      total_bytes += window_size_bytes(*input_grouped_block_residues_window_B);
    ASSERT(total_bytes <= max_shared_memory_bytes, DEBUG_STRING(total_bytes),
           DEBUG_STRING(max_shared_memory_bytes));
  }

// detect OpenBLAS using any OpenBLAS-specific macro from cblas.h
#ifdef OPENBLAS_THREAD
  // Disable BLAS threading explicitly, each rank should work single-threaded.
  // If this does not work, try also setting environment variables:
  // export OPENBLAS_NUM_THREADS=1
  // or also:
  // export OMP_NUM_THREADS=1
  openblas_set_num_threads(1);
#endif
  // TODO: disable threading for other CBLAS implementations too.
}

El::Int BigInt_Shared_Memory_Syrk_Context::input_group_height_per_prime() const
{
  return input_grouped_block_residues_window_A->block_residues.at(0)
    .at(group_index)
    .Height();
}
