#include "Block_Parser.hxx"
#include "../../../../SDP.hxx"
#include "../../../../../sdp_convert.hxx"

#include <boost/archive/binary_iarchive.hpp>
#include <rapidjson/istreamwrapper.h>

namespace
{
  void set_bilinear(const El::Matrix<El::BigFloat> &input,
                    El::Matrix<El::BigFloat> &local)
  {
    size_t height = input.Height();
    size_t width = input.Width() == 0 ? 1 : input.Width();
    local.Resize(height, width);
    for(size_t row(0); row != height; ++row)
      {
        for(size_t column(0); column != width; ++column)
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
    int width = value.at(0).size();
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

  // TODO move this code closer to Dual_Constraint_Group definition.
  Dual_Constraint_Group
  read_Dual_Constraint_Group(size_t index, std::istream &block_stream,
                             Block_File_Format format)
  {
    Dual_Constraint_Group group;
    if(format == bin)
      {
        // NB: this should match Dual_Constraint_Group.cxx - serialize_bin()
        boost::archive::binary_iarchive ar(block_stream);
        ar >> group;
      }
    else if(format == json)
      {
        rapidjson::IStreamWrapper wrapper(block_stream);
        Block_Parser parser;
        rapidjson::Reader reader;
        reader.Parse(wrapper, parser);

        group.block_index = index;
        group.dim = parser.dim;
        group.num_points = parser.num_points;
        group.constraint_matrix = to_matrix(parser.B_state);
        group.constraint_constants = parser.c_state.value;
        group.bilinear_bases[0] = to_matrix(parser.bilinear_bases_even_state);
        group.bilinear_bases[1] = to_matrix(parser.bilinear_bases_odd_state);
      }
    else
      {
        El::RuntimeError("Unknown Block_File_Format: ", format);
      }
    return group;
  }
}

void read_block_stream(
  const El::Grid &grid, const size_t &index, std::istream &block_stream,
  Block_File_Format format, SDP &sdp,
  std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local)
{
  auto group = read_Dual_Constraint_Group(index, block_stream, format);
  auto &c(sdp.primal_objective_c.blocks.at(index));
  c.SetGrid(grid);
  c.Resize(group.constraint_constants.size(), 1);
  size_t local_height(c.LocalHeight());
  if(c.GlobalCol(0) == 0)
    {
      for(size_t row = 0; row < local_height; ++row)
        {
          const size_t global_row(c.GlobalRow(row));
          c.SetLocal(row, 0, group.constraint_constants.at(global_row));
        }
    }

  set_bilinear(group.bilinear_bases[0], bilinear_bases_local.at(2 * index));
  set_bilinear(group.bilinear_bases[1],
               bilinear_bases_local.at(2 * index + 1));

  {
    auto &B(sdp.free_var_matrix.blocks.at(index));
    B.SetGrid(grid);
    B.Resize(group.constraint_matrix.Height(),
             group.constraint_matrix.Width());
    for(int64_t local_row(0); local_row != B.LocalHeight(); ++local_row)
      {
        const El::Int global_row(B.GlobalRow(local_row));
        for(int64_t local_column(0); local_column != B.LocalWidth();
            ++local_column)
          {
            const El::Int global_column(B.GlobalCol(local_column));
            B.SetLocal(local_row, local_column,
                       group.constraint_matrix.Get(global_row, global_column));
          }
      }
  }
}
