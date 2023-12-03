#include "sdp_convert/sdp_convert.hxx"
#include "sdp_read/sdp_read.hxx"

namespace fs = std::filesystem;

void read_xml_input(
  const fs::path &input_file, std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices,
  size_t &num_processed);

void read_pvm_input(
  const std::vector<fs::path> &input_files,
  std::vector<El::BigFloat> &dual_objectives_b,
  std::vector<Polynomial_Vector_Matrix> &polynomial_vector_matrices,
  size_t &num_processed)
{
  for(auto &input_file : input_files)
    {
      if(input_file.empty())
        {
          continue;
        }
      if(input_file.extension() == ".nsv")
        {
          read_pvm_input(read_nsv_file_list(input_file), dual_objectives_b,
                         polynomial_vector_matrices, num_processed);
        }
      else
        {
          read_xml_input(input_file, dual_objectives_b,
                         polynomial_vector_matrices, num_processed);
        }
    }
}
