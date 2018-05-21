#include "../Block_Diagonal_Matrix.hxx"

std::ostream &operator<<(std::ostream &os, const Block_Diagonal_Matrix &A)
{
  os << "{";
  for(auto block(A.blocks_elemental.begin());
      block != A.blocks_elemental.end();)
    {
      El::Print(*block,"",os);
      ++block;
      if(block != A.blocks_elemental.end())
        {
          os << ", ";
        }
    }
  os << "}";
  return os;
}
