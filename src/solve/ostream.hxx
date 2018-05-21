#pragma once

#include <ostream>
#include <vector>
#include <array>

template <class T>
std::ostream &operator<<(std::ostream &os, const std::array<T,2> &v)
{
  os << "{";
  for(auto element(v.begin()); element != v.end();)
    {
      os << *element;
      ++element;
      if(element != v.end())
        {
          os << ", ";
        }
    }
  os << "}";
  return os;
}

// print any vector<T>, including Vector
template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
{
  os << "{";
  for(auto element(v.begin()); element != v.end();)
    {
      os << *element;
      ++element;
      if(element != v.end())
        {
          os << ", ";
        }
    }
  os << "}";
  return os;
}
