#include "../../SDP_Solver.hxx"

// result = \sum_p a[p] A_p,
//
// where a[p] is a vector of length primalObjective.size() and the
// constraint matrices A_p are given by
//
//   A_(j,r,s,k) = \sum_{b \in blocks[j]}
//                     Block_b(v_{b,k} v_{b,k}^T \otimes E^{rs}),
//
// where v_{b,k} is the k-th column of bilinearBases[b], as described
// in SDP.h.
//
// Inputs: sdp, a
// Output: result (overwritten)
//

void constraint_matrix_weighted_sum(const SDP &sdp, const Block_Vector &a,
                                    Block_Diagonal_Matrix &result)
{
  for(size_t jj = 0; jj < sdp.dimensions.size(); ++jj)
    {
      const size_t block_size(sdp.degrees[jj] + 1);
      for(size_t block_index(2 * jj); block_index < 2 * jj + 2; ++block_index)
        {
          El::Zero(result.blocks[block_index]);
          for(size_t column_block = 0; column_block < sdp.dimensions[jj];
              ++column_block)
            for(size_t row_block = 0; row_block <= column_block; ++row_block)
              {
                const size_t result_block_size(
                  sdp.bilinear_bases_dist[block_index].Height());
                const size_t column_offset(column_block * result_block_size),
                  row_offset(row_block * result_block_size);
                size_t vector_offset(
                  ((column_block * (column_block + 1)) / 2 + row_block)
                  * block_size);
                El::DistMatrix<El::BigFloat> sub_vector(El::LockedView(
                  a.blocks[jj], vector_offset, 0, block_size, 1));
                El::DistMatrix<El::BigFloat> scaled_bases(
                  sdp.bilinear_bases_dist[block_index]);

                El::DiagonalScale(El::LeftOrRight::RIGHT,
                                  El::Orientation::NORMAL, sub_vector,
                                  scaled_bases);

                El::DistMatrix<El::BigFloat> result_sub_block(El::View(
                  result.blocks[block_index], row_offset, column_offset,
                  result_block_size, result_block_size));
                El::Gemm(El::Orientation::NORMAL, El::Orientation::TRANSPOSE,
                         El::BigFloat(column_block == row_block ? 1 : 0.5),
                         sdp.bilinear_bases_dist[block_index], scaled_bases,
                         El::BigFloat(0), result_sub_block);
              }
          if(sdp.dimensions[jj] > 1)
            {
              El::MakeSymmetric(El::UpperOrLowerNS::UPPER,
                                result.blocks[block_index]);
            }
        }
    }
}
