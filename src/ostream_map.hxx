#pragma once

#include "ostream_pair.hxx"
#include <map>

template <typename Key, typename Value>
std::ostream &operator<<(std::ostream &os, const std::map<Key, Value> &map)
{
  os << "[";
  auto iter(map.begin());
  if(iter!=map.end())
    {
      os << *iter;
      ++iter;
      for(; iter != map.end(); ++iter)
        {
          os << ", " << *iter;
        }
    }
  os << "]";
  return os;
}

