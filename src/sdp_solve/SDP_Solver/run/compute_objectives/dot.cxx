#include "sdp_solve/Block_Vector.hxx"
#include "sdpb_util/assert.hxx"

El::BigFloat dot(const Block_Vector &A, const Block_Vector &B)
{
  ASSERT(A.blocks.size() == B.blocks.size());
  El::BigFloat local_sum(0);
  for(size_t ii = 0; ii != A.blocks.size(); ++ii)
    {
      // FIXME: This feels slow.  It has to wait for each block
      // computation to be done before it can go to the next.
      local_sum += Dotu(A.blocks[ii], B.blocks[ii]);
    }
  // Make sure not to double count if blocks are distributed over more
  // than one processor.  We could also divide the sum by
  // X.blocks.front().Size().
  if(!A.blocks.empty() && A.blocks.front().Grid().Rank() != 0)
    {
      local_sum = 0;
    }
  return El::mpi::AllReduce(local_sum, El::mpi::COMM_WORLD);
}
