#pragma once

#include "block_timings.hxx"
#include "matrix_sizes.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_Comb.hxx"
#include "sdpb_util/Timers/Timers.hxx"

// Assign costs proportional to matrix sizes.
// This should balance out memory use when doing a timing run.
inline std::vector<Block_Cost>
block_costs_from_dimensions(const std::vector<size_t> &dimensions,
                            const std::vector<size_t> &num_points,
                            const size_t &dual_dimension)
{
  std::vector<Block_Cost> result;

  auto schur_sizes = schur_block_sizes(dimensions, num_points);
  const auto psd_sizes = psd_matrix_block_sizes(dimensions, num_points);
  const auto bilinear_sizes
    = bilinear_pairing_block_sizes(dimensions, num_points);

  auto elements_count
    = [](const std::vector<size_t> &sizes, const size_t index) {
        return sizes[index] * sizes[index];
      };

  // We store residues of L^{-1}B band modulo a bunch of primes
  // in a shared memory window of doubles,
  // see BigInt_Shared_Memory_Syrk_Context.input_block_residues_window.
  // residue_size is total size of residues for one BigFloat divided by size of BigFloat.
  double residue_size;
  {
    const auto total_schur_block_height
      = std::accumulate(schur_sizes.begin(), schur_sizes.end(), 0);

    // Same shared_memory_comm and Fmpz_Comb initialization,
    // as in initialize_bigint_syrk_context()
    El::mpi::Comm shared_memory_comm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                        &shared_memory_comm.comm);
    const auto num_nodes = El::mpi::Size() / El::mpi::Size(shared_memory_comm);
    const auto schur_block_height_per_node
      = total_schur_block_height / num_nodes;
    const Fmpz_Comb comb(El::gmp::Precision(), El::gmp::Precision(), 1,
                         schur_block_height_per_node);

    residue_size = static_cast<double>(comb.num_primes) * sizeof(double)
                   / El::BigFloat(1.1).SerializedSize();
    // Sanity check: residues always take more RAM
    // than the original BigFloat number,
    // otherwise the number cannot be restored.
    ASSERT(residue_size >= 1.0, residue_size);
  }

  for(size_t block = 0; block < schur_sizes.size(); ++block)
    {
      const auto schur = elements_count(schur_sizes, block);
      const auto psd = elements_count(psd_sizes, 2 * block)
                       + elements_count(psd_sizes, 2 * block + 1);
      const auto bilinear = elements_count(bilinear_sizes, 2 * block)
                            + elements_count(bilinear_sizes, 2 * block + 1);
      // P'xN, a band of B(=free_var_matrix)
      const auto B_band = schur_sizes[block] * dual_dimension;

      //  L^{-1}B residues, see BigInt_Shared_Memory_Syrk_Context.input_block_residues_window
      const size_t B_band_residues = std::round(B_band * residue_size);
      // Estimate total RAM associated with the block.
      // (There is also a RAM contribution from #(Q)=NxN, but it's
      // block-independent)
      auto total_cost
        = 2 * B_band + 5 * psd + 2 * schur + 2 * bilinear + B_band_residues;
      result.emplace_back(total_cost, block);
    }
  return result;
}
