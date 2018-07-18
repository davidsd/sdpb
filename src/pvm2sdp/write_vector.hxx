#pragma once

#include <boost/filesystem.hpp>
#include <vector>

template<typename T>
void write_vector(boost::filesystem::ofstream &output_stream,
                  const std::vector<T> &v)
{
  output_stream << v.size() << "\n";
  for(auto &element : v)
    {
      output_stream << element << "\n";
    }
}
