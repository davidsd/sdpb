#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"
#include "Block_File_Format.hxx"
#include "sdpb_util/assert.hxx"

#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <iostream>

namespace
{
  void write_matrix(std::ostream &output_stream,
                    const El::Matrix<El::BigFloat> &matrix,
                    const std::string &indentation)
  {
    output_stream << indentation << "[\n";
    for(int row = 0; row < matrix.Height(); ++row)
      {
        if(row != 0)
          {
            output_stream << ",\n";
          }
        output_stream << indentation << "  [\n";
        for(int column = 0; column < matrix.Width(); ++column)
          {
            if(column != 0)
              {
                output_stream << ",\n";
              }
            output_stream << indentation << "    \"" << matrix(row, column)
                          << "\"";
          }
        output_stream << "\n" << indentation << "  ]";
      }
    output_stream << "\n" << indentation << "]";
  }

  void write_bilinear_bases(std::ostream &output_stream,
                            const Dual_Constraint_Group &group)
  {
    output_stream << "  \"bilinear_bases_even\":\n";
    for(auto basis(group.bilinear_bases.begin());
        basis != group.bilinear_bases.end(); ++basis)
      {
        if(basis != group.bilinear_bases.begin())
          {
            output_stream << ",\n  \"bilinear_bases_odd\":\n";
          }
        // Ensure that each bilinearBasis is sampled the correct number
        // of times
        ASSERT(static_cast<size_t>(basis->Width()) == group.num_points);
        write_matrix(output_stream, *basis, "  ");
      }
    output_stream << ",\n";
  }

  void write_primal_objective_c(std::ostream &output_stream,
                                const Dual_Constraint_Group &group)
  {
    ASSERT(static_cast<size_t>(group.constraint_matrix.Height())
           == group.constraint_constants.size());

    output_stream << "  \"c\":\n";
    write_vector(output_stream, group.constraint_constants, "  ");
    output_stream << ",\n";
  }

  void write_free_var_matrix(std::ostream &output_stream,
                             const Dual_Constraint_Group &group)
  {
    output_stream << "  \"B\":\n";
    write_matrix(output_stream, group.constraint_matrix, "  ");
    output_stream << "\n";
  }

  void write_block_data_json(std::ostream &output_stream,
                             const Dual_Constraint_Group &group)
  {
    set_stream_precision(output_stream);
    output_stream << "{\n";
    write_bilinear_bases(output_stream, group);
    write_primal_objective_c(output_stream, group);
    write_free_var_matrix(output_stream, group);
    output_stream << "}\n";
  }

  void write_block_data_bin(std::ostream &output_stream,
                            const Dual_Constraint_Group &group)
  {
    boost::archive::binary_oarchive ar(output_stream);
    // store precision in order to ensure correct deserialization
    ar << El::gmp::Precision();
    // write fields in the same order as in Dual_Constraint_Group declaration
    // TODO use the same order in JSON?
    ar << group.constraint_matrix;
    ar << group.constraint_constants;
    ASSERT(group.bilinear_bases.size() == 2);
    ar << group.bilinear_bases[0];
    ar << group.bilinear_bases[1];
  }
}

void write_block_data(std::ostream &os, const Dual_Constraint_Group &group,
                      Block_File_Format format)
{
  switch(format)
    {
    case bin: write_block_data_bin(os, group); break;
    case json: write_block_data_json(os, group); break;
    default: RUNTIME_ERROR("Unknown Block_File_Format: ", format);
    }
}