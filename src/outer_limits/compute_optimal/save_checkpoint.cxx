#include "../../Verbosity.hxx"
#include "../../ostream_vector.hxx"
#include "../../ostream_set.hxx"
#include "../../set_stream_precision.hxx"

#include <El.hpp>

#include <boost/filesystem.hpp>
#include <boost/property_tree/json_parser.hpp>

void save_checkpoint(const boost::filesystem::path &checkpoint_directory,
                     const Verbosity &verbosity,
                     const boost::property_tree::ptree &parameter_properties,
                     const El::Matrix<El::BigFloat> &y,
                     const std::vector<std::set<El::BigFloat>> &points,
                     const El::BigFloat &infinity,
                     boost::optional<int64_t> &backup_generation,
                     int64_t &current_generation)
{
  if(!exists(checkpoint_directory))
    {
      create_directories(checkpoint_directory);
    }
  else if(!is_directory(checkpoint_directory))
    {
      throw std::runtime_error("Checkpoint directory '"
                               + checkpoint_directory.string()
                               + "'already exists, but is not a directory");
    }
  if(backup_generation)
    {
      remove(checkpoint_directory
             / ("points_" + std::to_string(backup_generation.value())
                + "_" + std::to_string(El::mpi::Rank()) + ".json"));
    }

  backup_generation = current_generation;
  current_generation += 1;
  boost::filesystem::path checkpoint_filename(
    checkpoint_directory
    / ("points_" + std::to_string(current_generation) + "_"
       + std::to_string(El::mpi::Rank()) + ".json"));

  const size_t max_retries(10);
  bool wrote_successfully(false);
  for(size_t attempt = 0; attempt < max_retries && !wrote_successfully;
      ++attempt)
    {
      boost::filesystem::ofstream checkpoint_stream(checkpoint_filename);
      set_stream_precision(checkpoint_stream);
      if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          std::cout << "Saving checkpoint to    : " << checkpoint_directory
                    << '\n';
        }
      checkpoint_stream << "{\n  \"points\":\n  [\n";
      // Output 'points' manually because we want to substitute in
      // 'inf' for infinity.
      for(auto block(points.begin()); block!=points.end(); ++block)
        {
          if(block != points.begin())
            {
              checkpoint_stream << ",\n";
            }
          checkpoint_stream << "    [\n";
          for(auto element(block->begin()); element!=block->end(); ++element)
            {
              if(element!=block->begin())
                {
                  checkpoint_stream << ",\n";
                }
              checkpoint_stream << "      \"";
              if(*element!=infinity)
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
      
      checkpoint_stream  << "\n  ]\n}\n";
      wrote_successfully = checkpoint_stream.good();
      if(!wrote_successfully)
        {
          if(attempt + 1 < max_retries)
            {
              std::stringstream ss;
              ss << "Error writing checkpoint file '" << checkpoint_filename
                 << "'.  Retrying " << (attempt + 2) << "/" << max_retries
                 << "\n";
              std::cerr << ss.str() << std::flush;
            }
          else
            {
              std::stringstream ss;
              ss << "Error writing checkpoint file '" << checkpoint_filename
                 << "'.  Exceeded max retries.\n";
              throw std::runtime_error(ss.str());
            }
        }
    }
  if(El::mpi::Rank() == 0)
    {
      boost::filesystem::ofstream y_stream(checkpoint_directory
                                           / "y_new.json");
      set_stream_precision(y_stream);
      y_stream << "{\n    \"current\": " << current_generation << ",\n"
               << "    \"backup\": " << backup_generation.value() << ",\n"
               << "    \"version\": \"" << SDPB_VERSION_STRING
               << "\",\n    \"options\": \n";

      boost::property_tree::write_json(y_stream, parameter_properties);
      y_stream << ",\n  \"y\":\n  [\n";
      for(int64_t row(0); row != y.Height(); ++row)
        {
          if(row != 0)
            {
              y_stream << ",\n";
            }
          y_stream << "    \"" << y(row, 0) << '"';
        }
      y_stream << "\n  ]\n}\n";
    }
  El::mpi::Barrier(El::mpi::COMM_WORLD);
  if(El::mpi::Rank() == 0)
    {
      rename(checkpoint_directory / "y_new.json",
             checkpoint_directory / "y.json");
    }
}
