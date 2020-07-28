#pragma once

#include <iostream>
#include <set>

template <template <typename> typename Sequence, typename T>
std::ostream &operator<<(std::ostream &os, const Sequence<T> &seq)
{
  os << "[";
  auto s(seq.begin());
  if(s!=seq.end())
    {
      os << *s;
      ++s;
      for(; s != seq.end(); ++s)
        {
          os << ", " << *s;
        }
    }
  os << "]";
  return os;
}

