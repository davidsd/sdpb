#pragma once

#include <iostream>
#include <utility>

template <typename T, typename U>
std::ostream &operator<<(std::ostream &os, const std::pair<T, U> &pair)
{
  os << "[" << pair.first << ", " << pair.second << "]";
  return os;
}

