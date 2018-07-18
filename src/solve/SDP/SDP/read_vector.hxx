#pragma once

#include <fstream>
#include <vector>

template<typename T>
void read_vector(std::ifstream &input_stream, std::vector<T> &v)
{
  size_t size;
  input_stream >> size;
  T element;
  for(size_t row=0; row<size; ++row)
    {
      input_stream >> element;
      v.push_back(element);
    }
}
