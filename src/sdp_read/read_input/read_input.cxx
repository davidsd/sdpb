#include "../Positive_Matrix_With_Prefactor.hxx"
#include "../../sdp_read.hxx"

#include <boost/filesystem.hpp>

void read_json(const boost::filesystem::path &input_path,
               std::vector<El::BigFloat> &objectives,
               std::vector<El::BigFloat> &normalization,
               std::vector<Positive_Matrix_With_Prefactor> &matrices);

void read_mathematica(const boost::filesystem::path &input_path,
                      std::vector<El::BigFloat> &objectives,
                      std::vector<El::BigFloat> &normalization,
                      std::vector<Positive_Matrix_With_Prefactor> &matrices);

void read_input(const boost::filesystem::path &input_file,
                std::vector<El::BigFloat> &objectives,
                std::vector<El::BigFloat> &normalization,
                std::vector<Positive_Matrix_With_Prefactor> &matrices)
{
  if(input_file.extension() == ".nsv")
    {
      for(auto &filename : read_file_list(input_file))
        {
          read_input(filename, objectives, normalization, matrices);
        }
    }
  else if(input_file.extension() == ".json")
    {
      read_json(input_file, objectives, normalization, matrices);
    }
  else
    {
      read_mathematica(input_file, objectives, normalization, matrices);
    }

  for(auto &matrix : matrices)
    {
      for(auto &pole : matrix.damped_rational.poles)
        {
          if(pole > 0)
            {
              throw std::runtime_error(
                "All poles must be positive, but found '" + to_string(pole)
                + "'");
            }
        }
    }
}
