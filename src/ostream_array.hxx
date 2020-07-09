#pragma once

#include <iostream>
#include <array>

template <typename T, size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T, N> &a)
{
  os << "[";
  auto v(a.begin());
  if(v!=a.end())
    {
      os << *v;
      ++v;
      for(; v != a.end(); ++v)
        {
          os << ", " << *v;
        }
    }
  os << "]";
  return os;
}

