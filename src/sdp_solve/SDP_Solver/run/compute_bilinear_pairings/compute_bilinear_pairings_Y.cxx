#include "../../../Block_Diagonal_Matrix.hxx"

// bilinear_pairings_Y[b] = Q[b]'^T A[b] Q[b]' for each block 0 <= b < Q.size()
// bilinear_pairings_Y[b], A[b] denote the b-th blocks of bilinear_pairings_Y,
// A, resp.

// Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product

// for each b, L.blocks[b], Q[b], Work[b], and bilinear_pairings_Y.blocks[b]
// must have the structure described above for `tensorTransposeCongruence'

void compute_bilinear_pairings_Y(
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  Block_Diagonal_Matrix &bilinear_pairings_Y)
{
  auto Y_block(Y.blocks.begin());
  auto bilinear_pairings_Y_block(bilinear_pairings_Y.blocks.begin());

  for(auto &block : bases_blocks)
    {
      auto temp_space(block);
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1),
           *Y_block, block, El::BigFloat(0), temp_space);
      Gemm(El::Orientation::TRANSPOSE, El::Orientation::NORMAL,
           El::BigFloat(1), block, temp_space, El::BigFloat(0),
           *bilinear_pairings_Y_block);
      El::MakeSymmetric(El::UpperOrLower::LOWER, *bilinear_pairings_Y_block);
      ++Y_block;
      ++bilinear_pairings_Y_block;
    }
}
