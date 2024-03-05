#include "SDP_Block_Data.hxx"

#include "Json_Block_Data_Parser.hxx"
#include "sdp_solve/SDP/SDP/set_bases_blocks.hxx"
#include "sdpb_util/Vector_State.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/boost_serialization.hxx"

#include <boost/archive/binary_iarchive.hpp>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/reader.h>

namespace
{
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
Block_Data_Parse_Result
parse_block_data(std::istream &block_stream, Block_File_Format format)
{
  Block_Data_Parse_Result result;
  if(format == bin)
    {
      // NB: this should match pmp2sdp/write_block_data.cxx
      boost::archive::binary_iarchive ar(block_stream);
      mp_bitcnt_t precision;
      ar >> precision;
      ASSERT(precision == El::gmp::Precision(),
             "Read GMP precision: ", precision,
             ", expected: ", El::gmp::Precision());
      ar >> result.B;
      ar >> result.c;
      ar >> result.bilinear_bases_even;
      ar >> result.bilinear_bases_odd;
    }
  else if(format == json)
    {
      rapidjson::IStreamWrapper wrapper(block_stream);
      Json_Block_Data_Parser parser(
        [&result](Block_Data_Parse_Result &&block_data_parse_result) {
          result = std::move(block_data_parse_result);
        });
      rapidjson::Reader reader;
      reader.Parse(wrapper, parser);
    }
  else
    {
      RUNTIME_ERROR("Unknown Block_File_Format: ", format);
    }
  return result;
}

SDP_Block_Data::SDP_Block_Data(std::istream &block_stream,
                               const Block_File_Format format,
                               const size_t block_index_local,
                               const Block_Info &block_info)
    : block_index_local(block_index_local)
{
  auto parse_result = parse_block_data(block_stream, format);
  constraint_matrix = std::move(parse_result.B);
  primal_objective_c = to_matrix(parse_result.c);
  bilinear_bases[0] = std::move(parse_result.bilinear_bases_even);
  bilinear_bases[1] = std::move(parse_result.bilinear_bases_odd);

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
