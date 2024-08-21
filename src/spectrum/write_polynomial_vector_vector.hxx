#pragma once

#include "pmp/Polynomial.hxx"
#include "pmp2sdp/write_vector.hxx"

void write_polynomial_vector_vector(std::ostream &output_stream,
                                    const std::vector<Polynomial_Vector> &m,
                                    const std::string &indentation)
{
  output_stream << indentation << "[\n";
  for(auto v(m.begin()); v != m.end(); ++v)
    {
      if(v != m.begin())
        {
          output_stream << ",\n";
        }

      // output_stream << indentation << "  \"" << *element << "\"";
      write_vector(output_stream,
                   static_cast<const std::vector<Polynomial> &>(*v),
                   indentation + "  ");
    }
  output_stream << "\n" << indentation << "]";
}