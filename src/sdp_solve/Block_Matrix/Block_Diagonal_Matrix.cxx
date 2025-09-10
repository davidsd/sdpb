#include "Block_Diagonal_Matrix.hxx"

Block_Diagonal_Matrix::Block_Diagonal_Matrix(
  const std::vector<size_t> &block_sizes,
  const std::vector<size_t> &block_indices, const El::Grid &grid)
{
  blocks.reserve(block_indices.size());
  for(auto &block_index : block_indices)
    {
      Abstract_Block_Diagonal_Matrix::add_block(block_sizes.at(block_index),
                                                grid);
    }
}

Paired_Block_Diagonal_Matrix::Paired_Block_Diagonal_Matrix(
  const std::vector<size_t> &block_sizes,
  const std::vector<size_t> &block_indices, const El::Grid &grid)
{
  blocks.reserve(block_indices.size() * 2);
  for(auto &block_index : block_indices)
    {
      Abstract_Block_Diagonal_Matrix::add_block(
        block_sizes.at(block_index * 2), grid);
      Abstract_Block_Diagonal_Matrix::add_block(
        block_sizes.at(block_index * 2 + 1), grid);
    }
}
El::DistMatrix<El::BigFloat> &
Paired_Block_Diagonal_Matrix::get_block(const size_t block_index,
                                        const size_t parity)
{
  ASSERT(parity == 0 || parity == 1, DEBUG_STRING(parity));
  return blocks.at(2 * block_index + parity);
}
const El::DistMatrix<El::BigFloat> &
Paired_Block_Diagonal_Matrix::get_block(const size_t block_index,
                                        const size_t parity) const
{
  ASSERT(parity == 0 || parity == 1, DEBUG_STRING(parity));
  return blocks.at(2 * block_index + parity);
}
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L,
                            const std::vector<size_t> &block_indices,
                            const std::string &name)
{
  ASSERT_EQUAL(A.blocks.size(), L.blocks.size(), name);
  ASSERT_EQUAL(A.blocks.size(), block_indices.size(), name);
  for(size_t block = 0; block < block_indices.size(); ++block)
    {
      const size_t block_index = block_indices[block];
      auto &L_block = L.blocks.at(block);
      auto &A_block = A.blocks.at(block);

      // FIXME: Use pivoting?
      L_block = A_block;
      try
        {
          Cholesky(El::UpperOrLowerNS::LOWER, L_block);
        }
      catch(std::exception &e)
        {
          RUNTIME_ERROR("Error when computing Cholesky decomposition of "
                        "Paired_Block_Diagonal_Matrix ",
                        name, DEBUG_STRING(block_index), ": ", e.what());
        }
    }
}
void cholesky_decomposition(const Paired_Block_Diagonal_Matrix &A,
                            Paired_Block_Diagonal_Matrix &L,
                            const std::vector<size_t> &block_indices,
                            const std::string &name)
{
  ASSERT_EQUAL(A.blocks.size(), L.blocks.size(), name);
  ASSERT_EQUAL(A.blocks.size(), 2 * block_indices.size(), name);
  for(size_t block = 0; block < block_indices.size(); ++block)
    {
      const size_t block_index = block_indices[block];
      for(const size_t parity : {0, 1})
        {
          auto &L_block = L.get_block(block, parity);
          auto &A_block = A.get_block(block, parity);

          // FIXME: Use pivoting?
          L_block = A_block;
          try
            {
              Cholesky(El::UpperOrLowerNS::LOWER, L_block);
            }
          catch(std::exception &e)
            {
              RUNTIME_ERROR("Error when computing Cholesky decomposition of "
                            "Paired_Block_Diagonal_Matrix ",
                            name, DEBUG_STRING(block_index),
                            DEBUG_STRING(parity), ": ", e.what());
            }
        }
    }
}
