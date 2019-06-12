#include "../../SDP_Solver.hxx"

// dualResidues[p] = primalObjective[p] - Tr(A_p Y) - (FreeVarMatrix y)_p,
// for 0 <= p < primalObjective.size()
//
// The pairings Tr(A_p Y) can be written in terms of A_Y:
//
//   Tr(A_(j,r,s,k) Y) = \sum_{b \in blocks[j]}
//                       (1/2) (A_Y_{ej r + k, ej s + k} +
//                              swap (r <-> s))
// where ej = d_j + 1.

void compute_dual_residues_and_error(const Block_Info &block_info,
                                     const SDP &sdp, const Block_Vector &y,
                                     const Block_Diagonal_Matrix &A_Y,
                                     Block_Vector &dual_residues,
                                     El::BigFloat &dual_error, Timers &timers)
{
  auto &dual_residues_timer(timers.add_and_start("run.computeDualResidues"));

  auto dual_residues_block(dual_residues.blocks.begin());
  auto primal_objective_c_block(sdp.primal_objective_c.blocks.begin());
  auto y_block(y.blocks.begin());
  auto B_block(sdp.B.blocks.begin());
  auto A_Y_block(A_Y.blocks.begin());

  El::BigFloat local_max(0);
  for(auto &block_index : block_info.block_indices)
    {
      Zero(*dual_residues_block);
      const size_t block_size(block_info.degrees[block_index] + 1);

      // Not sure whether it is better to first loop over blocks in
      // the result or over sub-blocks in A_Y
      for(size_t bb = 2 * block_index; bb < 2 * block_index + 2; ++bb)
        {
          for(size_t column_block = 0;
              column_block < block_info.dimensions[block_index];
              ++column_block)
            for(size_t row_block = 0; row_block <= column_block; ++row_block)
              {
                size_t column_offset(column_block * block_size),
                  row_offset(row_block * block_size);

                El::DistMatrix<El::BigFloat> lower_diagonal(El::GetDiagonal(
                  El::LockedView(*A_Y_block, row_offset, column_offset,
                                 block_size, block_size)));

                size_t residue_row_offset(
                  ((column_block * (column_block + 1)) / 2 + row_block)
                  * block_size);

                El::DistMatrix<El::BigFloat> residue_sub_block(El::View(
                  *dual_residues_block, residue_row_offset, 0, block_size, 1));

                // Do a little less work when on a diagonal block.
                if(column_offset == row_offset)
                  {
                    El::Axpy(El::BigFloat(-1.0), lower_diagonal,
                             residue_sub_block);
                  }
                else
                  {
                    El::Axpy(El::BigFloat(-0.5), lower_diagonal,
                             residue_sub_block);

                    El::DistMatrix<El::BigFloat> upper_diagonal(
                      El::GetDiagonal(El::LockedView(*A_Y_block, column_offset,
                                                     row_offset, block_size,
                                                     block_size)));

                    El::Axpy(El::BigFloat(-0.5), upper_diagonal,
                             residue_sub_block);
                  }
              }
          ++A_Y_block;
        }
      // dualResidues -= FreeVarMatrix * y
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(-1),
           *B_block, *y_block, El::BigFloat(1), *dual_residues_block);
      // dualResidues += primalObjective
      Axpy(El::BigFloat(1), *primal_objective_c_block, *dual_residues_block);

      local_max = Max(local_max, El::MaxAbs(*dual_residues_block));

      ++primal_objective_c_block;
      ++y_block;
      ++B_block;
      ++dual_residues_block;
    }
  dual_error
    = El::mpi::AllReduce(local_max, El::mpi::MAX, El::mpi::COMM_WORLD);
  dual_residues_timer.stop();
}
