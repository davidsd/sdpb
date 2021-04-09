#pragma once

#include <iostream>
#include <vector>

template <typename T>
inline void write_vector(std::ostream &output_stream, const std::vector<T> &v,
                         const std::string &indentation)
{
  output_stream << indentation << "[\n";
  for(auto element(v.begin()); element != v.end(); ++element)
    {
      if(element != v.begin())
        {
          output_stream << ",\n";
        }
      output_stream << indentation << "  \"" << *element << "\"";
    }
  output_stream << "\n" << indentation << "]";
}

template <typename T>
inline void
write_vector(std::ostream &output_stream, const std::vector<T> &v,
             const std::string &indentation, const std::string &name)
{
  output_stream << indentation << "\"" << name << "\":\n";
  write_vector(output_stream, v, indentation);
}

inline void
write_vector(std::ostream &output_stream, const std::vector<size_t> &v)
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

inline void write_vector(std::ostream &output_stream,
                         const std::vector<size_t> &v, const std::string &name)
{
  output_stream << "\"" << name << "\": ";
  write_vector(output_stream, v);
}
