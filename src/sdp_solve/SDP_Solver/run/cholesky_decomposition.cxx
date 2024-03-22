#include "sdp_solve/Block_Diagonal_Matrix.hxx"
#include "sdp_solve/Block_Info.hxx"

// Compute L (lower triangular) such that A = L L^T
void cholesky_decomposition(const Block_Diagonal_Matrix &A,
                            Block_Diagonal_Matrix &L,
                            const Block_Info &block_info,
                            const std::string &name)
{
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      // FIXME: Use pivoting?
      L.blocks[b] = A.blocks[b];
      try
        {
          Cholesky(El::UpperOrLowerNS::LOWER, L.blocks[b]);
        }
      catch(std::exception &e)
        {
          RUNTIME_ERROR("Error when computing Cholesky decomposition of "
                        "Block_Diagonal_Matrix ",
                        name,
                        ", block index = ", block_info.block_indices.at(b),
                        ": ", e.what());
        }
    }
}
