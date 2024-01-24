#include "SDP_Block_Data.hxx"

#include "Block_Parser/Block_Parser.hxx"
#include "sdp_solve/SDP/SDP/set_bases_blocks.hxx"
#include "sdpb_util/Number_State.hxx"
#include "sdpb_util/Vector_State.hxx"
#include "sdpb_util/boost_serialization.hxx"

#include <boost/archive/binary_iarchive.hpp>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/reader.h>

namespace
{
  template <class FloatType>
  El::Matrix<FloatType>
  to_matrix(const Vector_State<Vector_State<Number_State<FloatType>>> &input)
  {
    const auto &value = input.value;
    int height = value.size();
    int width = value.empty() ? 1 : value.at(0).size();
    El::Matrix<FloatType> result(height, width);

    for(int row = 0; row < height; ++row)
      {
        for(int col = 0; col < width; ++col)
          {
            result.Set(row, col, value.at(row).at(col));
          }
      }
    return result;
  }

  // Convert vec[N] to Matrix(N,1)
  template <class FloatType>
  El::Matrix<FloatType> to_matrix(const std::vector<FloatType> &input)
  {
    int height = input.size();
    int width = 1;
    El::Matrix<FloatType> result(height, width);
    for(int i = 0; i < height; ++i)
      {
        result.Set(i, 0, input.at(i));
      }
    return result;
  }
}

// TODO move this code closer to Dual_Constraint_Group definition?
void parse_block_data(std::istream &block_stream, Block_File_Format format,
                      El::Matrix<El::BigFloat> &constraint_matrix,
                      std::vector<El::BigFloat> &constraint_constants,
                      El::Matrix<El::BigFloat> &bilinear_bases_even,
                      El::Matrix<El::BigFloat> &bilinear_bases_odd)
{
  if(format == bin)
    {
      // NB: this should match pmp2sdp/write_block_data.cxx
      boost::archive::binary_iarchive ar(block_stream);
      mp_bitcnt_t precision;
      ar >> precision;
      if(precision != El::gmp::Precision())
        {
          El::RuntimeError("Read GMP precision: ", precision,
                           ", expected: ", El::gmp::Precision());
        }
      ar >> constraint_matrix;
      ar >> constraint_constants;
      ar >> bilinear_bases_even;
      ar >> bilinear_bases_odd;
    }
  else if(format == json)
    {
      rapidjson::IStreamWrapper wrapper(block_stream);
      Block_Parser parser;
      rapidjson::Reader reader;
      reader.Parse(wrapper, parser);

      constraint_matrix = to_matrix(parser.B_state);
      constraint_constants = parser.c_state.value;
      bilinear_bases_even = to_matrix(parser.bilinear_bases_even_state);
      bilinear_bases_odd = to_matrix(parser.bilinear_bases_odd_state);
    }
  else
    {
      El::RuntimeError("Unknown Block_File_Format: ", format);
    }
}

SDP_Block_Data::SDP_Block_Data(std::istream &block_stream,
                               const Block_File_Format format,
                               const size_t block_index_local,
                               const Block_Info &block_info)
    : block_index_local(block_index_local)
{
  std::vector<El::BigFloat> constraint_constants;
  parse_block_data(block_stream, format, constraint_matrix,
                   constraint_constants, bilinear_bases[0], bilinear_bases[1]);
  primal_objective_c = to_matrix(constraint_constants);

  const size_t block_index = block_info.block_indices.at(block_index_local);

  // Set bases_blocks:
  for(const size_t parity : {0, 1})
    {
      const auto &bilinear_base_local = bilinear_bases[parity];
      auto &bases_block = bases_blocks[parity];

      const auto height
        = block_info.get_psd_matrix_block_size(block_index, parity);
      const auto width
        = block_info.get_bilinear_pairing_block_size(block_index, parity);
      bases_block.Resize(height, width);

      set_bilinear_bases_block_local(bilinear_base_local, bases_block);
    }
}
