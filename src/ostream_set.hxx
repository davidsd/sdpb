#pragma once

#include <iostream>
#include <set>

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::set<T> &vec)
{
  os << "[";
  auto v(vec.begin());
  if(v!=vec.end())
    {
      os << *v;
      ++v;
      for(; v != vec.end(); ++v)
        {
          os << ", " << *v;
        }
    }
  os << "]";
  return os;
}

