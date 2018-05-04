#include "../../SDP_Solver.hxx"

// dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
// for 0 <= p < primalObjective.size()
//
// The pairings Tr(A_p Y) can be written in terms of BilinearPairingsY:
//
//   Tr(A_(j,r,s,k) Y) = \sum_{b \in blocks[j]}
//                       (1/2) (BilinearPairingsY_{ej r + k, ej s + k} +
//                              swap (r <-> s))
// where ej = d_j + 1.
//
// Inputs: sdp, y, BilinearPairingsY
// Output: dualResidues (overwriten)
//
void compute_dual_residues(const SDP &sdp, const Vector &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Vector &dual_residues)
{
  for(size_t j = 0; j < sdp.dimensions.size(); j++)
    {
      const int ej = sdp.degrees[j] + 1;

      for(auto &t : sdp.constraint_indices[j])
        {
          const int p = t.p;
          const int ej_r = t.r * ej;
          const int ej_s = t.s * ej;
          const int k = t.k;

          // dualResidues[p] = -Tr(A_p Y)
          dual_residues[p] = 0;
          for(auto &b : sdp.blocks[j])
            {
              dual_residues[p]
                -= bilinear_pairings_Y.blocks[b].elt(ej_r + k, ej_s + k);
              dual_residues[p]
                -= bilinear_pairings_Y.blocks[b].elt(ej_s + k, ej_r + k);
            }
          dual_residues[p] /= 2;

          // dualResidues[p] = -Tr(A_p Y) - (FreeVarMatrix y)_p
          for(size_t n = 0; n < sdp.free_var_matrix.cols; n++)
            {
              dual_residues[p] -= sdp.free_var_matrix.elt(p, n) * y[n];
            }

          // dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix
          // y)_p
          dual_residues[p] += sdp.primal_objective_c[p];
        }
    }
}

void compute_dual_residues(const SDP &sdp,
                           const El::DistMatrix<El::BigFloat> &y,
                           const Block_Diagonal_Matrix &bilinear_pairings_Y,
                           Block_Matrix &dual_residues)
{
  // Maps from linear subblock offsets in dual_residues to subblock
  // indices (r,s) for bilinear_pairings_Y.
  //
  // This is only valid up to 5x5 blocks, but it should fail noisily if it is
  // out of bounds.
  const std::vector<size_t> r_map(
    {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4}),
    s_map({0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4});

  // Set up all of the queues to fetch the diagonals and off-diagonals
  // of bilinear_pairings_Y for Tr(A_p Y).  For small blocks, this
  // is probably all local.
  for(size_t jj = 0; jj < sdp.dimensions.size(); ++jj)
    {
      const size_t bilinear_bases_block_index(2 * jj);
      const size_t bilinear_bases_row_size(
        sdp.bilinear_bases_elemental_local[bilinear_bases_block_index].Height()
        * sdp.dimensions[jj]);

      const size_t row_offset(dual_residues.blocks[jj].GlobalRow(0));
      const size_t local_height(dual_residues.blocks[jj].LocalHeight());
      for(size_t row = 0; row < local_height; ++row)
        {
          const size_t global_row(row + row_offset);
          const size_t subblock_index(global_row / bilinear_bases_row_size),
            k(global_row % bilinear_bases_row_size);
          const size_t r(r_map.at(subblock_index)),
            s(s_map.at(subblock_index));
          const size_t r_k(r * bilinear_bases_row_size + k),
            s_k(s * bilinear_bases_row_size + k);

          bilinear_pairings_Y.blocks_elemental[bilinear_bases_block_index]
            .QueuePull(r_k, s_k);
          bilinear_pairings_Y.blocks_elemental[bilinear_bases_block_index]
            .QueuePull(s_k, r_k);
          bilinear_pairings_Y.blocks_elemental[bilinear_bases_block_index + 1]
            .QueuePull(r_k, s_k);
          bilinear_pairings_Y.blocks_elemental[bilinear_bases_block_index + 1]
            .QueuePull(s_k, r_k);
        }
    }

  // dualResidues = primalObjective - (FreeVarMatrix y)
  for(size_t b = 0; b < dual_residues.blocks.size(); ++b)
    {
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(-1),
           sdp.free_var_matrix_elemental.blocks[b], y, El::BigFloat(0),
           dual_residues.blocks[b]);
      Axpy(El::BigFloat(1), sdp.primal_objective_c_elemental.blocks[b],
           dual_residues.blocks[b]);
    }

  // dualResidues -= Tr(A_p Y)
  //
  // Doing this later helps overlap communication and computation
  // (maybe?)
  std::vector<El::BigFloat> even_queue, odd_queue;
  for(size_t jj = 0; jj < sdp.dimensions.size(); ++jj)
    {
      const size_t bilinear_bases_block_index(2 * jj);
      bilinear_pairings_Y.blocks_elemental[bilinear_bases_block_index]
        .ProcessPullQueue(even_queue);
      bilinear_pairings_Y.blocks_elemental[bilinear_bases_block_index+1]
        .ProcessPullQueue(odd_queue);

      const size_t local_height(dual_residues.blocks[jj].LocalHeight());
      for(size_t row = 0; row < local_height; ++row)
        {
          const size_t queue_index(2 * row);
          dual_residues.blocks[jj].UpdateLocal(
            row, 0,
            (even_queue.at(queue_index) + even_queue.at(queue_index + 1)
             + odd_queue.at(queue_index) + odd_queue.at(queue_index + 1))
              / -2);
        }
    }
}
