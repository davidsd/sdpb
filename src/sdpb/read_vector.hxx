#pragma once

#include <fstream>
#include <vector>

template <typename T>
void read_vector(std::ifstream &input_stream, std::vector<T> &v)
{
  size_t size;
  input_stream >> size;
  if(!input_stream.good())
    {
      throw std::runtime_error("Error reading vector size from file");
    }
  v.reserve(size);
  T element;
  for(size_t row = 0; row < size; ++row)
    {
      input_stream >> element;
      v.push_back(element);
    }
  if(!input_stream.good())
    {
      throw std::runtime_error("Error reading vector elements from file");
    }
}
