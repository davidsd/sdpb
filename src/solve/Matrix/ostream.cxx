#include "../Matrix.hxx"

std::ostream &operator<<(std::ostream &os, const Matrix &a)
{
  os << "{";
  for(int r = 0; r < a.rows; r++)
    {
      os << "{";
      for(int c = 0; c < a.cols; c++)
        {
          os << a.elt(r, c);
          if(c < a.cols - 1)
            {
              os << ", ";
            }
        }
      os << "}";
      if(r < a.rows - 1)
        {
          os << ", ";
        }
    }
  os << "}";
  return os;
}
