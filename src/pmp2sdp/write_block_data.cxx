#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <cmath>
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
        ASSERT_EQUAL(basis->Width(), group.num_points);
        write_matrix(output_stream, *basis, "  ");
      }
    output_stream << ",\n";
  }

  void write_primal_objective_c(std::ostream &output_stream,
                                const Dual_Constraint_Group &group)
  {
    ASSERT_EQUAL(group.constraint_matrix.Height(),
                 group.constraint_constants.size());

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

  // pmp2sdp --encodingPoints computes at a higher precision than it writes
  // (set_block_data_output_precision); 0: the working precision.
  mp_bitcnt_t block_data_output_precision = 0;
  mp_bitcnt_t output_precision()
  {
    return block_data_output_precision != 0 ? block_data_output_precision
                                            : El::gmp::Precision();
  }

  void write_block_data_json(std::ostream &output_stream,
                             const Dual_Constraint_Group &group)
  {
    set_stream_precision(output_stream);
    if(block_data_output_precision != 0)
      output_stream.precision(
        std::ceil(block_data_output_precision * std::log10(2.0)) + 1);
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
    // (the numbers themselves must have been truncated to this precision,
    // see set_block_data_output_precision)
    ar << output_precision();
    // write fields in the same order as in Dual_Constraint_Group declaration
    // TODO use the same order in JSON?
    ar << group.constraint_matrix;
    ar << group.constraint_constants;
    ASSERT_EQUAL(group.bilinear_bases.size(), 2);
    ar << group.bilinear_bases[0];
    ar << group.bilinear_bases[1];
  }
}

void set_block_data_output_precision(const mp_bitcnt_t precision)
{
  // GMP rounds a precision up to whole limbs (224 -> 256): the block files
  // carry the rounded value, as El::gmp::Precision() would report it.
  block_data_output_precision = 0;
  if(precision != 0)
    {
      mpf_t rounded;
      mpf_init2(rounded, precision);
      block_data_output_precision = mpf_get_prec(rounded);
      mpf_clear(rounded);
    }
  // Binary blocks: El's layout with the output precision's limb count
  // (sdpb_util/boost_serialization.hxx)
  sdpb_serialization::output_num_limbs()
    = block_data_output_precision == 0
          || block_data_output_precision == El::gmp::Precision()
        ? 0
        : sdpb_serialization::num_limbs_for_precision(
          block_data_output_precision);
}

void write_block_data(std::ostream &os, const Dual_Constraint_Group &group,
                      Block_File_Format format)
{
  switch(format)
    {
    case Block_File_Format::bin: write_block_data_bin(os, group); break;
    case Block_File_Format::json: write_block_data_json(os, group); break;
    default: RUNTIME_ERROR("Unknown Block_File_Format: ", format);
    }
}