#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"
#include "sdpb_util/json/Json_Writer.hxx"
#include "spectrum/Zeros.hxx"

#include <fstream>

namespace fs = std::filesystem;

void write_file(const fs::path &output_path,
                const std::vector<Zeros> &zeros_blocks, Timers &timers)
{
  if(El::mpi::Rank() != 0)
    return;

  Scoped_Timer timer(timers, "write_file");

  if(output_path.has_parent_path())
    fs::create_directories(output_path.parent_path());
  std::ofstream ofs(output_path);
  ASSERT(ofs.good(), "Problem when opening output file: ", output_path);

  rapidjson::OStreamWrapper ows(ofs);
  Json_PrettyWriter writer(ows);
  writer.SetIndent(' ', 2);

  writer.StartArray();
  for(auto zeros_iterator(zeros_blocks.begin());
      zeros_iterator != zeros_blocks.end(); ++zeros_iterator)
    {
      const auto &zeros = *zeros_iterator;
      const auto block_path = zeros_iterator->block_path.string();
      // TODO store block_index in Zeros?
      ASSERT(!block_path.empty(), "Empty path for block_",
             std::distance(zeros_blocks.begin(), zeros_iterator));
      writer.StartObject();
      {
        writer.Key("block_path");
        writer.String(block_path.c_str());

        writer.Key("zeros");
        writer.StartArray();
        for(const auto &zero : zeros.zeros)
          {
            writer.StartObject();
            {
              writer.Key("zero");
              writer.BigFloat(zero.zero);

              writer.Key("lambda");
              const auto &lambda = zero.lambda;
              writer.StartArray();
              // Print each column, i.e. each eigenvector
              for(int col = 0; col < lambda.Width(); ++col)
                {
                  writer.StartArray();
                  for(int row = 0; row < lambda.Height(); ++row)
                    writer.BigFloat(lambda(row, col));
                  writer.EndArray();
                }
              writer.EndArray();
            }
            writer.EndObject();
          }
        writer.EndArray();

        writer.Key("error");
        writer.BigFloat(zeros.error);
      }
      writer.EndObject();
    }
  writer.EndArray();
  ASSERT(writer.IsComplete());
  ASSERT(ofs.good(), "Problem when writing to output file: ", output_path);
}
