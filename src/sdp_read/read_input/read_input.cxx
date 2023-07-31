#include "../Positive_Matrix_With_Prefactor.hxx"
#include "../../sdp_read.hxx"

#include <boost/filesystem.hpp>

void read_json(const boost::filesystem::path &input_path,
               std::vector<El::BigFloat> &objectives,
               std::vector<El::BigFloat> &normalization,
               std::vector<Positive_Matrix_With_Prefactor> &matrices,
               size_t &num_processed);

void read_mathematica(const boost::filesystem::path &input_path,
                      std::vector<El::BigFloat> &objectives,
                      std::vector<El::BigFloat> &normalization,
                      std::vector<Positive_Matrix_With_Prefactor> &matrices,
                      size_t &num_processed);

void read_input(const boost::filesystem::path &input_file,
                std::vector<El::BigFloat> &objectives,
                std::vector<El::BigFloat> &normalization,
                std::vector<Positive_Matrix_With_Prefactor> &matrices,
                size_t &num_processed)
{
  if(!boost::filesystem::exists(input_file))
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
  else
    {
      read_mathematica(input_file, objectives, normalization, matrices,
                       num_processed);
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
