#include "Checkpoint_Parser.hxx"

#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"
#include "sdpb_util/ostream/ostream_set.hxx"

#include <El.hpp>

#include <rapidjson/istreamwrapper.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/optional.hpp>
#include <filesystem>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <algorithm>

namespace fs = std::filesystem;

void load_checkpoint(
  const fs::path &checkpoint_directory, const Verbosity &verbosity, boost::optional<int64_t> &backup_generation,
  int64_t &current_generation,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &yp_to_y_star,
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &dual_objective_b_star,
  El::Matrix<El::BigFloat> &yp, std::vector<std::set<El::BigFloat>> &points,
  El::BigFloat &threshold, El::BigFloat &primal_c_scale)
{
  if(fs::exists(checkpoint_directory)
     && fs::is_directory(checkpoint_directory))
    {
      fs::directory_iterator entry(checkpoint_directory);
      int64_t max_generation(-1);
      const std::string prefix("checkpoint_"), suffix(".json.gz");
      for(fs::directory_iterator entry(checkpoint_directory);
          entry != fs::directory_iterator(); ++entry)
        {
          const std::string name(entry->path().filename().string());
          if(boost::algorithm::starts_with(name, "checkpoint_"))
            {
              max_generation
                = std::max(max_generation,
                           int64_t(std::stoll(name.substr(prefix.size()))));
            }
        }
      if(max_generation != -1)
        {
          backup_generation = max_generation;
          current_generation = max_generation + 1;

          const std::string checkpoint_name(
            (checkpoint_directory
             / (prefix + std::to_string(max_generation) + suffix))
              .string());

          if(El::mpi::Rank() == 0 && verbosity >= Verbosity::regular)
            {
              std::cout << "Loading checkpoint: " << checkpoint_name << "\n";
            }

          boost::iostreams::filtering_istream input_stream;
          input_stream.push(boost::iostreams::gzip_decompressor());
          input_stream.push(boost::iostreams::file_source(checkpoint_name));

          rapidjson::IStreamWrapper wrapper(input_stream);
          Checkpoint_Parser parser;
          rapidjson::Reader reader;
          reader.Parse(wrapper, parser);

          ASSERT(parser.generation_state.value == El::BigFloat(max_generation),
                 "Invalid checkpoint. The generation from the file name ",
                 max_generation,
                 " does not match the generation stored inside: ",
                 parser.generation_state.value);

          yp_to_y_star.Resize(parser.y_transform_state.value.size(),
                              parser.y_transform_state.value.size());
          for(int64_t row(0); row < yp_to_y_star.LocalHeight(); ++row)
            {
              for(int64_t column(0); column < yp_to_y_star.LocalWidth();
                  ++column)
                {
                  yp_to_y_star.SetLocal(
                    row, column,
                    parser.y_transform_state.value.at(row).at(column));
                }
            }

          dual_objective_b_star.Resize(parser.b_state.value.size(), 1);
          for(int64_t row(0); row < dual_objective_b_star.Height(); ++row)
            {
              dual_objective_b_star.SetLocal(row, 0,
                                             parser.b_state.value[row]);
            }

          yp.Resize(parser.yp_state.value.size(), 1);
          for(int64_t row(0); row < yp.Height(); ++row)
            {
              yp(row, 0) = parser.yp_state.value[row];
            }

          points.clear();
          for(auto &block : parser.points_state.value)
            {
              points.emplace_back();
              auto &back(points.back());
              for(auto &point : block)
                {
                  back.emplace(point);
                }
            }
          threshold = parser.threshold_state.value;
          primal_c_scale = parser.c_scale_state.value;
        }
    }
}
