#include "../Block_Diagonal_Matrix.hxx"

// Minimum eigenvalue of A, via the QR method
// Inputs:
// A           : symmetric Block_Diagonal_Matrix
//
El::BigFloat min_eigenvalue(Block_Diagonal_Matrix &A)
{
  El::BigFloat lambda_min(El::limits::Max<El::BigFloat>());

  // FIXME: El::HermitianEig() seems broken in parallel.  So this
  // makes a local copy of the blocks before computing the
  // eigenvalues.  So it is reeeeeeeealy inefficient.

  // El::DistMatrix<El::BigFloat> eigenvalues;
  El::Matrix<El::BigFloat> local_eigenvalues;
  for(auto &block : A.blocks_elemental)
    {
      El::Matrix<El::BigFloat> local_copy(block.Height(), block.Width());
      block.ReservePulls(local_copy.Height() * local_copy.Width());

      for(size_t row = 0; row < local_copy.Height(); ++row)
        for(size_t column = 0; column < local_copy.Width(); ++column)
          {
            block.QueuePull(row, column);
          }
      std::vector<El::BigFloat> pull_queue;
      block.ProcessPullQueue(pull_queue);

      size_t queue_index(0);
      for(size_t row = 0; row < local_copy.Height(); ++row)
        for(size_t column = 0; column < local_copy.Width(); ++column)
          {
            local_copy.Set(row, column, pull_queue[queue_index]);
            ++queue_index;
          }
      El::HermitianEig(El::UpperOrLowerNS::LOWER, local_copy,
                       local_eigenvalues);
      lambda_min = El::Min(lambda_min, El::Min(local_eigenvalues));

      /// FIXME: Maybe use regular double precision here?
      // El::HermitianEig(El::UpperOrLowerNS::LOWER, block, eigenvalues);
      // lambda_min = El::Min(lambda_min, El::Min(eigenvalues));
    }
  return lambda_min;
}
