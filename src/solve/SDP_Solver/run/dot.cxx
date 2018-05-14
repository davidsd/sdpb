#include "../../Block_Vector.hxx"
#include <cassert>

El::BigFloat dot(const Block_Vector &a, const Block_Vector &b)
{
  assert(a.blocks.size() == b.blocks.size());
  El::BigFloat result(0);
  for(size_t ii=0; ii!=a.blocks.size(); ++ii)
    {
      // FIXME: This feels slow.  It has to wait for each block
      // computation to be done before it can go to the next.
      result+=Dotu(a.blocks[ii],b.blocks[ii]);
    }
  return result;
}
