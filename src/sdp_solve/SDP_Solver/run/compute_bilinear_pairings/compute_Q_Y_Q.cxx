#include "../../../Block_Diagonal_Matrix.hxx"

// Q_Y_Q[b] = Q[b]'^T A[b] Q[b]' for each block 0 <= b < Q.size()
// Q_Y_Q[b], A[b] denote the b-th blocks of Q_Y_Q,
// A, resp.

// Q[b]' = Q[b] \otimes 1, where \otimes denotes tensor product

// for each b, L.blocks[b], Q[b], Work[b], and Q_Y_Q.blocks[b]
// must have the structure described above for `tensorTransposeCongruence'

void compute_Q_Y_Q(
  const Block_Diagonal_Matrix &Y,
  const std::vector<El::DistMatrix<El::BigFloat>> &bases_blocks,
  Block_Diagonal_Matrix &Q_Y_Q)
{
  auto Y_block(Y.blocks.begin());
  auto Q_Y_Q_block(Q_Y_Q.blocks.begin());

  for(auto &block : bases_blocks)
    {
      auto temp_space(block);
      Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, El::BigFloat(1),
           *Y_block, block, El::BigFloat(0), temp_space);
      Gemm(El::Orientation::TRANSPOSE, El::Orientation::NORMAL,
           El::BigFloat(1), block, temp_space, El::BigFloat(0),
           *Q_Y_Q_block);
      El::MakeSymmetric(El::UpperOrLower::LOWER, *Q_Y_Q_block);
      ++Y_block;
      ++Q_Y_Q_block;
    }
}
