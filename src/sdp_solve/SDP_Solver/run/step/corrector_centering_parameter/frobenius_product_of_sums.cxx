#include "sdp_solve/Block_Diagonal_Matrix.hxx"

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
//
El::BigFloat frobenius_product_of_sums(const Block_Diagonal_Matrix &X,
                                       const Block_Diagonal_Matrix &dX,
                                       const Block_Diagonal_Matrix &Y,
                                       const Block_Diagonal_Matrix &dY)
{
  El::BigFloat local_sum(0);
  for(size_t b = 0; b < X.blocks.size(); b++)
    {
      // FIXME: This can be sped up by not have intermediate results.
      // It may require looking into the implementation of Dotu.
      El::DistMatrix<El::BigFloat> X_dX(X.blocks[b]);
      X_dX += dX.blocks[b];
      El::DistMatrix<El::BigFloat> Y_dY(Y.blocks[b]);
      Y_dY += dY.blocks[b];
      local_sum += Dotu(X_dX, Y_dY);
    }

  // Make sure not to double count if blocks are distributed over more
  // than one processor.  We could also divide the sum by
  // X.blocks.front().Size().
  if(!X.blocks.empty() && X.blocks.front().Grid().Rank() != 0)
    {
      local_sum = 0;
    }
  return El::mpi::AllReduce(local_sum, El::mpi::COMM_WORLD);
}
