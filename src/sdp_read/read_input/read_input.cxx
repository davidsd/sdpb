#include "sdp_read/sdp_read.hxx"
#include "sdp_read/Positive_Matrix_With_Prefactor.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void read_json(const fs::path &input_path,
               std::vector<El::BigFloat> &objectives,
               std::vector<El::BigFloat> &normalization,
               std::vector<Positive_Matrix_With_Prefactor> &matrices,
               size_t &num_processed);

void read_mathematica(const fs::path &input_path,
                      std::vector<El::BigFloat> &objectives,
                      std::vector<El::BigFloat> &normalization,
                      std::vector<Positive_Matrix_With_Prefactor> &matrices,
                      size_t &num_processed);

void read_input(const fs::path &input_file,
                std::vector<El::BigFloat> &objectives,
                std::vector<El::BigFloat> &normalization,
                std::vector<Positive_Matrix_With_Prefactor> &matrices,
                size_t &num_processed)
{
  if(!exists(input_file))
    {
      El::RuntimeError("Cannot find input file: ", input_file);
    }
  if(input_file.extension() == ".nsv")
    {
      for(auto &filename : read_nsv_file_list(input_file))
        {
          read_input(filename, objectives, normalization, matrices,
                     num_processed);
        }
    }
  else if(input_file.extension() == ".json")
    {
      read_json(input_file, objectives, normalization, matrices,
                num_processed);
    }
  else if(input_file.extension() == ".m")
    {
      read_mathematica(input_file, objectives, normalization, matrices,
                       num_processed);
    }
  else
    {
      El::RuntimeError("Cannot parse input file: ", input_file,
                       ". Expected .nsv, .json or .m extension.");
    }

  for(auto &matrix : matrices)
    {
      for(auto &pole : matrix.damped_rational.poles)
        {
          if(pole > 0)
            {
              throw std::runtime_error(
                "All poles must be negative, but found '" + to_string(pole)
                + "'");
            }
        }
    }
}
