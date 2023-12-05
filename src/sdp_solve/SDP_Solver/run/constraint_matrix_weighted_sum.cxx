#include "sdp_solve/SDP_Solver.hxx"

// result = \sum_p a[p] A_p,
//
// where a[p] is a vector of length primalObjective.size() and the
// constraint matrices A_p are given by
//
//   A_(j,r,s,k) = \sum_{b \in blocks[j]}
//                     Block_b(v_{b,k} v_{b,k}^T \otimes E^{rs}),
//
// where v_{b,k} is the k-th column of bilinear_bases[b], as described
// in SDP.h.

void constraint_matrix_weighted_sum(const Block_Info &block_info,
                                    const SDP &sdp, const Block_Vector &a,
                                    Block_Diagonal_Matrix &result)
{
  auto a_block(a.blocks.begin());
  auto result_block(result.blocks.begin());
  auto bilinear_bases_block(sdp.bilinear_bases.begin());

  for(auto &block_index : block_info.block_indices)
    {
      const size_t block_size(block_info.num_points[block_index]);
      for(size_t parity = 0; parity < 2; ++parity)
        {
          // TODO: Remove this because it gets zero'd out in Gemm?
          El::Zero(*result_block);
          for(size_t column_block = 0;
              column_block < block_info.dimensions[block_index];
              ++column_block)
            for(size_t row_block = 0; row_block <= column_block; ++row_block)
              {
                const size_t result_block_size(bilinear_bases_block->Height());
                const size_t column_offset(column_block * result_block_size),
                  row_offset(row_block * result_block_size);
                size_t vector_offset(
                  ((column_block * (column_block + 1)) / 2 + row_block)
                  * block_size);
                El::DistMatrix<El::BigFloat> sub_vector(
                  El::LockedView(*a_block, vector_offset, 0, block_size, 1));
                El::DistMatrix<El::BigFloat> scaled_bases(
                  *bilinear_bases_block);

                El::DiagonalScale(El::LeftOrRight::RIGHT,
                                  El::Orientation::NORMAL, sub_vector,
                                  scaled_bases);

                El::DistMatrix<El::BigFloat> result_sub_block(
                  El::View(*result_block, row_offset, column_offset,
                           result_block_size, result_block_size));
                El::Gemm(El::Orientation::NORMAL, El::Orientation::TRANSPOSE,
                         El::BigFloat(column_block == row_block ? 1 : 0.5),
                         *bilinear_bases_block, scaled_bases, El::BigFloat(0),
                         result_sub_block);
              }
          if(block_info.dimensions[block_index] > 1)
            {
              El::MakeSymmetric(El::UpperOrLowerNS::UPPER, *result_block);
            }
          ++result_block;
          ++bilinear_bases_block;
        }
      ++a_block;
    }
}
