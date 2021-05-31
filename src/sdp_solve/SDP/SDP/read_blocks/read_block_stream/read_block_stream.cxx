#include "Block_Parser.hxx"
#include "../../../../SDP.hxx"

#include <rapidjson/istreamwrapper.h>

namespace
{
  void set_bilinear(const std::vector<std::vector<El::BigFloat>> &input,
                    El::Matrix<El::BigFloat> &local)
  {
    const size_t height(input.size()),
      width(input.empty() ? 1 : input.at(0).size());
    local.Resize(height, width);
    for(size_t row(0); row != height; ++row)
      {
        for(size_t column(0); column != width; ++column)
          {
            local(row, column) = input[row].at(column);
          }
      }
  }
}

void read_block_stream(
  const El::Grid &grid, const size_t &index, std::istream &block_stream,
  SDP &sdp, std::vector<El::Matrix<El::BigFloat>> &bilinear_bases_local)
{
  rapidjson::IStreamWrapper wrapper(block_stream);
  Block_Parser parser;
  rapidjson::Reader reader;
  reader.Parse(wrapper, parser);

  auto &c(sdp.primal_objective_c.blocks.at(index));
  c.SetGrid(grid);
  c.Resize(parser.c_state.value.size(), 1);
  size_t local_height(c.LocalHeight());
  if(c.GlobalCol(0) == 0)
    {
      for(size_t row = 0; row < local_height; ++row)
        {
          const size_t global_row(c.GlobalRow(row));
          c.SetLocal(row, 0, parser.c_state.value.at(global_row));
        }
    }

  set_bilinear(parser.bilinear_bases_even_state.value,
               bilinear_bases_local.at(2 * index));
  set_bilinear(parser.bilinear_bases_odd_state.value,
               bilinear_bases_local.at(2 * index + 1));

  {
    auto &B(sdp.free_var_matrix.blocks.at(index));
    B.SetGrid(grid);
    B.Resize(parser.B_state.value.size(), parser.B_state.value.at(0).size());
    for(int64_t local_row(0); local_row != B.LocalHeight(); ++local_row)
      {
        const El::Int global_row(B.GlobalRow(local_row));
        for(int64_t local_column(0); local_column != B.LocalWidth();
            ++local_column)
          {
            const El::Int global_column(B.GlobalCol(local_column));
            B.SetLocal(local_row, local_column,
                       parser.B_state.value.at(global_row).at(global_column));
          }
      }
  }
}
