#include "../BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"
#include "sdpb_util/ostream/pretty_print_bytes.hxx"

#include <El.hpp>

#include <cblas.h>

BigInt_Shared_Memory_Syrk_Context::BigInt_Shared_Memory_Syrk_Context(
  const Bigint_Syrk_Config &cfg, const size_t group_index,
  const std::vector<size_t> &block_index_local_to_global,
  const Verbosity verbosity,
  const std::function<Job_Schedule(Blas_Job::Kind kind, El::UpperOrLower uplo,
                                   size_t num_ranks, size_t num_primes,
                                   int output_height, int output_width,
                                   Verbosity _verbosity)> &create_job_schedule)
    : shared_memory_comm(cfg.shared_memory_comm),
      group_index(group_index),
      num_groups(cfg.input_height_per_group.size()),
      total_block_height_per_node(cfg.input_height),
      comb(cfg.comb()),
      verbosity(verbosity),
      input_window_split_factor(cfg.input_split_factor),
      output_window_split_factor(cfg.output_split_factor),
      block_index_local_to_global(block_index_local_to_global),
      create_blas_job_schedule_func(create_job_schedule)
{
  size_t window_width = cfg.window_width();
  input_grouped_block_residues_window_A
    = std::make_unique<Vertical_Block_Matrix_Residues_Window<double>>(
      shared_memory_comm, comb.num_primes, num_groups,
      cfg.input_window_height_per_group_per_prime, window_width);
  if(cfg.num_input_windows() != 1)
    {
      ASSERT_EQUAL(cfg.num_input_windows(), 2);
      input_grouped_block_residues_window_B
        = std::make_unique<Vertical_Block_Matrix_Residues_Window<double>>(
          shared_memory_comm, comb.num_primes, num_groups,
          cfg.input_window_height_per_group_per_prime, window_width);
    }
  output_residues_window = std::make_unique<Matrix_Residues_Window<double>>(
    shared_memory_comm, comb.num_primes, window_width, window_width);
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

  // Sanity check
  ASSERT_EQUAL(cfg.input_window_size(),
               input_grouped_block_residues_window_A->size());
  ASSERT_EQUAL(cfg.output_window_size(), output_residues_window->size());

  // Print sizes
  if(verbosity >= Verbosity::debug && shared_memory_comm.Rank() == 0)
    {
      std::ostringstream os;
      El::BuildStream(os, "create BigInt_Shared_Memory_Syrk_Context, rank=",
                      El::mpi::Rank(), "\n");
      El::BuildStream(os, "  Number of primes: ", comb.num_primes, "\n");
      El::BuildStream(os, "  Blocks on the node:\n");
      El::BuildStream(os, "    Total height: ",
                      static_cast<const size_t>(total_block_height_per_node),
                      "\n");
      El::BuildStream(os, "    Width: ", cfg.input_width, "\n");
      El::BuildStream(os, "    Elements: ",
                      total_block_height_per_node * cfg.input_width, "\n");
      // TODO for some reason, El::BuildStream() fails for std::vector
      // although we've included sdpb_util/ostream/ostream_vector.hxx
      os << "    Heights for each MPI group: " << cfg.input_height_per_group
         << "\n";

      El::BuildStream(os, "  Output residues window (Q):\n");
      El::BuildStream(os, "    Window size: ",
                      pretty_print_bytes(
                        output_residues_window->size() * sizeof(double), true),
                      "\n");
      El::BuildStream(os, "    Split factor: ", output_window_split_factor,
                      "\n");
      El::BuildStream(os, "    Height=Width (per prime): ", window_width,
                      "\n");
      El::BuildStream(os, "    Total elements (per prime): ",
                      window_width * window_width, "\n");

      auto input_window_height = cfg.input_window_height();

      El::BuildStream(os, "  Input residues window (P):\n");
      El::BuildStream(os, "    Number of windows: ", cfg.num_input_windows(),
                      "\n");
      El::BuildStream(
        os, "    Window size: ",
        pretty_print_bytes(input_grouped_block_residues_window_A->size()
                             * sizeof(double),
                           true),
        "\n");
      El::BuildStream(os, "    Split factor: ", input_window_split_factor,
                      "\n");
      El::BuildStream(os, "    Height (per prime): ", input_window_height,
                      "\n");
      El::BuildStream(os, "    Width: ", window_width, "\n");
      os << "    Heights for each MPI group: "
         << cfg.input_window_height_per_group_per_prime << "\n";

      if(cfg.get_reduce_scatter_buffer_size() > 0)
        {
          El::BuildStream(
            os, "  Buffer size for reduce-scatter Q: ",
            pretty_print_bytes(cfg.get_reduce_scatter_buffer_size(), true),
            "\n");
        }

      El::Output(os.str());
    }
}

El::Int BigInt_Shared_Memory_Syrk_Context::input_group_height_per_prime() const
{
  return input_grouped_block_residues_window_A->block_residues.at(0)
    .at(group_index)
    .Height();
}
