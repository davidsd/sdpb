#include "Block_Parser.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdp_convert/sdp_convert.hxx"

#include <boost/archive/binary_iarchive.hpp>
#include <rapidjson/istreamwrapper.h>

namespace
{
  void set_bilinear(const El::Matrix<El::BigFloat> &input,
                    El::Matrix<El::BigFloat> &local)
  {
    El::Int height = input.Height();
    El::Int width = input.Width();
    local.Resize(height, width);
    for(El::Int row(0); row != height; ++row)
      {
        for(El::Int column(0); column != width; ++column)
          {
            local(row, column) = input.Get(row, column);
          }
      }
  }

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

void read_block_stream(
  const El::Grid &grid, const size_t &index, std::istream &block_stream,
  Block_File_Format format, SDP &sdp,
  std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local)
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

  set_bilinear(bilinear_bases_even, bilinear_bases_local.at(2 * index));
  set_bilinear(bilinear_bases_odd, bilinear_bases_local.at(2 * index + 1));

  {
    auto &B(sdp.free_var_matrix.blocks.at(index));
    B.SetGrid(grid);
    B.Resize(constraint_matrix.Height(), constraint_matrix.Width());
    for(El::Int local_row(0); local_row != B.LocalHeight(); ++local_row)
      {
        const El::Int global_row(B.GlobalRow(local_row));
        for(El::Int local_column(0); local_column != B.LocalWidth();
            ++local_column)
          {
            const El::Int global_column(B.GlobalCol(local_column));
            B.SetLocal(local_row, local_column,
                       constraint_matrix.Get(global_row, global_column));
          }
      }
  }
}
