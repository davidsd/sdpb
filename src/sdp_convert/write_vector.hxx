#pragma once

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <vector>

template<typename T>
void write_vector(boost::filesystem::ofstream &output_stream,
                  const std::vector<T> &v, const std::string &name)
{
  output_stream << "\"" << name << "\": [";
  for(auto element(v.begin()); element!=v.end(); ++element)
    {
      if(element!=v.begin())
        {
          output_stream << ", ";
        }
      output_stream << "\"" << *element << "\"";
    }
  output_stream << "]";
}
