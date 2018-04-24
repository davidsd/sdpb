#include "../Matrix.hxx"

std::ostream &operator<<(std::ostream &os, const Matrix &a)
{
  os << "{";
  for(size_t r = 0; r < a.rows; r++)
    {
      os << "{";
      for(size_t c = 0; c < a.cols; c++)
        {
          os << a.elt(r, c);
          if(c + 1 < a.cols)
            {
              os << ", ";
            }
        }
      os << "}";
      if(r + 1 < a.rows)
        {
          os << ", ";
        }
    }
  os << "}";
  return os;
}
