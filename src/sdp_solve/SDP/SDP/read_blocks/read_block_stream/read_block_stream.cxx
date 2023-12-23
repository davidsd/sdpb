#include "Block_Parser.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdp_convert/sdp_convert.hxx"
#include "sdp_solve/SDP/SDP/set_bases_blocks.hxx"
#include "sdpb_util/copy_matrix.hxx"

#include <boost/archive/binary_iarchive.hpp>
#include <rapidjson/istreamwrapper.h>

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
      // NB: this should match sdp_convert/write_block_data.cxx
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

void read_block_stream(const El::Grid &grid, const size_t &index,
                       std::istream &block_stream, Block_File_Format format,
                       const Block_Info &block_info, SDP &sdp)
{
  El::Matrix<El::BigFloat> constraint_matrix;
  std::vector<El::BigFloat> constraint_constants;
  El::Matrix<El::BigFloat> bilinear_bases_even;
  El::Matrix<El::BigFloat> bilinear_bases_odd;

  parse_block_data(block_stream, format, constraint_matrix,
                   constraint_constants, bilinear_bases_even,
                   bilinear_bases_odd);

  auto &c(sdp.primal_objective_c.blocks.at(index));
  c.SetGrid(grid);
  c.Resize(constraint_constants.size(), 1);
  El::Int local_height(c.LocalHeight());
  if(c.GlobalCol(0) == 0)
    {
      for(El::Int row = 0; row < local_height; ++row)
        {
          const size_t global_row(c.GlobalRow(row));
          c.SetLocal(row, 0, constraint_constants.at(global_row));
        }
    }

  // sdp.bilinear_bases and sdp.bases_blocks
  for(const size_t parity : {0, 1})
    {
      const auto &bilinear_bases_local
        = parity == 0 ? bilinear_bases_even : bilinear_bases_odd;

      // Set sdp.bilinear_bases:

      const size_t bilinear_index_local = 2 * index + parity;
      auto &bilinear_base = sdp.bilinear_bases.at(bilinear_index_local);
      bilinear_base.SetGrid(grid);

      copy_matrix(bilinear_bases_local, bilinear_base);

      // Set sdp.bases_blocks:

      auto &bases_block = sdp.bases_blocks.at(bilinear_index_local);
      auto pairing_sizes(block_info.bilinear_pairing_block_sizes());
      auto psd_sizes(block_info.psd_matrix_block_sizes());

      size_t block_index = block_info.block_indices.at(index);
      El::Int height
        = block_info.get_psd_matrix_block_size(block_index, parity);
      El::Int width
        = block_info.get_bilinear_pairing_block_size(block_index, parity);
      bases_block.SetGrid(grid);
      bases_block.Resize(height, width);

      set_bilinear_bases_block(bilinear_bases_local, bases_block);
    }

  {
    auto &B(sdp.free_var_matrix.blocks.at(index));
    B.SetGrid(grid);
    B.Resize(constraint_matrix.Height(), constraint_matrix.Width());
    copy_matrix(constraint_matrix, B);
  }
}
