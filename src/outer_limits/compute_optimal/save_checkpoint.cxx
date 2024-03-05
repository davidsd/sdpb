#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"
#include "sdpb_util/ostream/ostream_set.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>

#include <boost/optional.hpp>
#include <filesystem>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace fs = std::filesystem;

namespace
{
  inline fs::path checkpoint_name(const int64_t &generation)
  {
    return "checkpoint_" + std::to_string(generation) + ".json.gz";
  }
}

void save_checkpoint(
  const fs::path &checkpoint_directory, const Verbosity &verbosity,
  const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &yp_to_y_star,
  const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &dual_objective_b_star,
  const El::Matrix<El::BigFloat> &yp,
  const std::vector<std::set<El::BigFloat>> &points,
  const El::BigFloat &infinity, const El::BigFloat &threshold,
  const El::BigFloat &primal_c_scale,
  boost::optional<int64_t> &backup_generation, int64_t &current_generation)
{
  if(checkpoint_directory.empty())
    {
      return;
    }
  if(El::mpi::Rank() == 0)
    {
      if(!exists(checkpoint_directory))
        {
          create_directories(checkpoint_directory);
        }
      else if(!is_directory(checkpoint_directory))
        {
          RUNTIME_ERROR("Checkpoint directory ", checkpoint_directory,
                        " already exists, but is not a directory");
        }
      if(backup_generation)
        {
          remove(checkpoint_directory
                 / checkpoint_name(backup_generation.value()));
        }

      backup_generation = current_generation;
      current_generation += 1;
      fs::path checkpoint_filename(checkpoint_directory / checkpoint_name(current_generation));

      const size_t max_retries(10);
      bool wrote_successfully(false);
      for(size_t attempt = 0; attempt < max_retries && !wrote_successfully;
          ++attempt)
        {
          boost::iostreams::filtering_ostream checkpoint_stream;
          // Use no compression.  Filter is for CRC only.
          checkpoint_stream.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(0)));
          checkpoint_stream.push(
            boost::iostreams::file_sink(checkpoint_filename.string()));

          set_stream_precision(checkpoint_stream);
          if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
            {
              std::cout << "Saving checkpoint to    : " << checkpoint_directory
                        << '\n';
            }

          checkpoint_stream << "{\n  \"generation\": \"" << current_generation
                            << "\",\n"
                            << "  \"threshold\": \"" << threshold << "\",\n"
                            << "  \"c_scale\": \"" << primal_c_scale << "\",\n"
                            << "  \"yp\":\n  [\n";
          for(int64_t row(0); row != yp.Height(); ++row)
            {
              if(row != 0)
                {
                  checkpoint_stream << ",\n";
                }
              checkpoint_stream << "    \"" << yp(row, 0) << '"';
            }
          checkpoint_stream << "\n  ],\n";
          checkpoint_stream << "  \"points\":\n  [\n";
          // Output 'points' manually because we want to substitute in
          // 'inf' for infinity.
          for(auto block(points.begin()); block != points.end(); ++block)
            {
              if(block != points.begin())
                {
                  checkpoint_stream << ",\n";
                }
              checkpoint_stream << "    [\n";
              for(auto element(block->begin()); element != block->end();
                  ++element)
                {
                  if(element != block->begin())
                    {
                      checkpoint_stream << ",\n";
                    }
                  checkpoint_stream << "      \"";
                  if(*element != infinity)
                    {
                      checkpoint_stream << *element;
                    }
                  else
                    {
                      checkpoint_stream << "inf";
                    }
                  checkpoint_stream << '"';
                }
              checkpoint_stream << "\n    ]";
            }
          checkpoint_stream << "\n  ],\n";

          checkpoint_stream << "  \"y_transform\":\n  [\n";
          for(int64_t row(0); row < yp_to_y_star.LocalHeight(); ++row)
            {
              if(row != 0)
                {
                  checkpoint_stream << ",\n";
                }
              checkpoint_stream << "    [\n";
              for(int64_t column(0); column < yp_to_y_star.LocalWidth();
                  ++column)
                {
                  if(column != 0)
                    {
                      checkpoint_stream << ",\n";
                    }
                  checkpoint_stream
                    << "      \"" << yp_to_y_star.GetLocal(row, column) << '"';
                }
              checkpoint_stream << "\n    ]";
            }
          checkpoint_stream << "\n  ],\n";

          checkpoint_stream << "  \"b\":\n  [\n";
          for(int64_t row(0); row < dual_objective_b_star.LocalHeight(); ++row)
            {
              if(row != 0)
                {
                  checkpoint_stream << ",\n";
                }
              checkpoint_stream
                << "      \"" << dual_objective_b_star.GetLocal(row, 0) << '"';
            }
          checkpoint_stream << "\n  ]\n}\n";

          wrote_successfully = checkpoint_stream.good();
          if(!wrote_successfully)
            {
              if(attempt + 1 < max_retries)
                {
                  std::stringstream ss;
                  ss << "Error writing checkpoint file '"
                     << checkpoint_filename << "'.  Retrying " << (attempt + 2)
                     << "/" << max_retries << "\n";
                  std::cerr << ss.str() << std::flush;
                }
              else
                {
                  RUNTIME_ERROR("Error writing checkpoint file ",
                                checkpoint_filename,
                                "'.  Exceeded max retries.");
                }
            }
        }
    }
}
