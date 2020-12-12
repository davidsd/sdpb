#pragma once

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <vector>


template <typename T>
inline void write_vector(boost::filesystem::ofstream &output_stream,
                         const std::vector<T> &v)
{
  output_stream << "[";
  for(auto element(v.begin()); element != v.end(); ++element)
    {
      if(element != v.begin())
        {
          output_stream << ", ";
        }
      output_stream << "\"" << *element << "\"";
    }
  output_stream << "]";
}

template <>
inline void write_vector(boost::filesystem::ofstream &output_stream,
                         const std::vector<size_t> &v)
{
  output_stream << "[";
  for(auto element(v.begin()); element != v.end(); ++element)
    {
      if(element != v.begin())
        {
          output_stream << ", ";
        }
      output_stream << *element;
    }
  output_stream << "]";
}

template <typename T>
inline void write_vector(boost::filesystem::ofstream &output_stream,
                         const std::vector<T> &v, const std::string &name)
{
  output_stream << "\"" << name << "\": ";
  write_vector(output_stream,v);
}
