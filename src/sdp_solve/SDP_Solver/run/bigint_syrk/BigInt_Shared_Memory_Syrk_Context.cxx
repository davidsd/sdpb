#include "BigInt_Shared_Memory_Syrk_Context.hxx"
#include "blas_jobs/Blas_Job_Schedule.hxx"
#include "fmpz/fmpz_mul_blas_util.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <cblas.h>

// code adopted from flint mul_blas.c
namespace
{
  El::Int sum(const std::vector<El::Int> &block_heights)
  {
    return std::accumulate(block_heights.begin(), block_heights.end(), 0);
  }

  template <class T>
  size_t window_size_bytes(const Residue_Matrices_Window<T> &window)
  {
    return window.height * window.width * window.num_primes * sizeof(T);
  }
}

BigInt_Shared_Memory_Syrk_Context::BigInt_Shared_Memory_Syrk_Context(
  const El::mpi::Comm &shared_memory_comm, const size_t group_index,
  const std::vector<int> &group_comm_sizes, const mp_bitcnt_t precision,
  const size_t max_shared_memory_bytes,
  const std::vector<El::Int> &blocks_height_per_group, const int block_width,
  const std::vector<size_t> &block_index_local_to_global, const bool debug,
  const std::function<Blas_Job_Schedule(size_t num_ranks, size_t num_primes,
                                        int output_width, bool debug)>
    &create_job_schedule)
    : shared_memory_comm(shared_memory_comm),
      group_index(group_index),
      group_comm_sizes(group_comm_sizes),
      num_groups(group_comm_sizes.size()),
      total_block_height_per_node(sum(blocks_height_per_group)),
      comb(precision, precision, 1, total_block_height_per_node),
      output_residues_window(shared_memory_comm, comb.num_primes, block_width,
                             block_width, debug),
      block_index_local_to_global(block_index_local_to_global),
      blas_job_schedule(create_job_schedule(
        shared_memory_comm.Size(), comb.num_primes, block_width, debug))
{
  ASSERT_EQUAL(blocks_height_per_group.size(), num_groups);

  // TODO split output window too
  const size_t output_window_bytes = window_size_bytes(output_residues_window);
  ASSERT(output_window_bytes < max_shared_memory_bytes,
         DEBUG_STRING(max_shared_memory_bytes),
         DEBUG_STRING(output_window_bytes));

  const size_t max_input_window_bytes
    = max_shared_memory_bytes - output_window_bytes;

  const size_t residues_bytes_per_single_block_row
    = block_width * comb.num_primes * sizeof(double);
  const size_t total_input_residues_bytes
    = total_block_height_per_node * residues_bytes_per_single_block_row;

  std::vector<int> input_window_height_per_group_per_prime(num_groups);
  if(total_input_residues_bytes <= max_input_window_bytes)
    {
      // no need to split
      input_window_height_per_group_per_prime = blocks_height_per_group;
      input_window_split_factor = 1;
    }
  else
    {
      const size_t input_window_height_per_prime
        = max_input_window_bytes / residues_bytes_per_single_block_row;
      for(size_t i = 0; i < num_groups; ++i)
        {
          input_window_height_per_group_per_prime.at(i)
            = input_window_height_per_prime * group_comm_sizes.at(i)
              / shared_memory_comm.Size();
          ASSERT(
            input_window_height_per_group_per_prime.at(i) > 0,
            "Cannot allocated shared memory window for block residues."
            " Required: at least ",
            residues_bytes_per_single_block_row,
            " bytes per each MPI group, available: ", max_input_window_bytes,
            " bytes in total. ", DEBUG_STRING(max_shared_memory_bytes),
            DEBUG_STRING(output_window_bytes));

          size_t curr_split_factor
            = blocks_height_per_group.at(i)
              / input_window_height_per_group_per_prime.at(i);
          if(blocks_height_per_group.at(i)
               % input_window_height_per_group_per_prime.at(i)
             != 0)
            {
              ++curr_split_factor;
            }
          // Set the same split factor for all ranks on a node,
          // so that bigint_syrk_blas_shmem() has the same number of iterations
          // (otherwise MPI will fail)
          input_window_split_factor
            = std::max(input_window_split_factor, curr_split_factor);
        }
    }
  input_grouped_block_residues_window
    = std::make_unique<Block_Residue_Matrices_Window<double>>(
      shared_memory_comm, comb.num_primes, num_groups,
      input_window_height_per_group_per_prime, block_width, debug);

  ASSERT(input_window_split_factor > 0);
  ASSERT(window_size_bytes(*input_grouped_block_residues_window)
           <= max_input_window_bytes,
         DEBUG_STRING(window_size_bytes(*input_grouped_block_residues_window)),
         DEBUG_STRING(max_input_window_bytes));

  // Disable BLAS threading explicitly, each rank should work single-threaded
  openblas_set_num_threads(1);

  // Clear input and outpur residue windows
  clear_residues();

  // Print sizes
  if(debug && shared_memory_comm.Rank() == 0)
    {
      // TODO print window sizes and split factors
      std::ostringstream os;
      El::BuildStream(os, "create BigInt_Shared_Memory_Syrk_Context, rank=",
                      El::mpi::Rank(), "\n");
      El::BuildStream(
        os, "  Shared memory limit, bytes: ", max_shared_memory_bytes, "\n");
      El::BuildStream(os, "  Number of primes: ", comb.num_primes, "\n");

      El::BuildStream(os, "  Blocks on the node:\n");
      El::BuildStream(os, "    Total height: ",
                      static_cast<const size_t>(total_block_height_per_node),
                      "\n");
      El::BuildStream(os, "    Width: ", block_width, "\n");
      El::BuildStream(
        os, "    Elements: ", total_block_height_per_node * block_width, "\n");

      El::BuildStream(os, "  Output residues window (Q):\n");
      El::BuildStream(os, "    Window size, bytes: ",
                      window_size_bytes(output_residues_window), "\n");
      El::BuildStream(os, "    Total elements (per prime): ",
                      block_width * block_width, "\n");

      El::BuildStream(os, "  Input residues window (P):\n");
      El::BuildStream(os, "    Window size, bytes: ",
                      window_size_bytes(*input_grouped_block_residues_window),
                      "\n");
      El::BuildStream(os, "    Split factor: ", input_window_split_factor,
                      "\n");
      El::BuildStream(os, "    Height (per prime): ",
                      input_grouped_block_residues_window->height, "\n");
      El::Print(input_window_height_per_group_per_prime,
                "    Heights for each MPI group:", ", ", os);
      os << "\n";
      El::Print(group_comm_sizes, "    MPI group sizes:", ", ", os);
      os << "\n";

      El::Output(os.str());
    }
}

El::Int BigInt_Shared_Memory_Syrk_Context::input_group_height_per_prime() const
{
  return input_grouped_block_residues_window->block_residues.at(0)
    .at(group_index)
    .Height();
}
