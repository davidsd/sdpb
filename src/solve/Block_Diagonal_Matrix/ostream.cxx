#include "../Block_Diagonal_Matrix.hxx"

std::ostream &operator<<(std::ostream &os, const Block_Diagonal_Matrix &A)
{
  os << "{";
  for(auto block(A.blocks.begin()); block != A.blocks.end();)
    {
      El::Print(*block, "", os);
      ++block;
      if(block != A.blocks.end())
        {
          os << ", ";
        }
    }
  os << "}";
  return os;
}
