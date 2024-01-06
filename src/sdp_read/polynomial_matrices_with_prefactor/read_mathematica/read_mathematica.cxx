#include "read_mathematica.hxx"

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

namespace fs = std::filesystem;

const char *
parse_SDP(const char *begin, const char *end,
          const std::function<bool(size_t matrix_index)> &should_parse_matrix,
          std::vector<El::BigFloat> &objectives,
          std::vector<El::BigFloat> &normalization, size_t &num_matrices,
          std::map<size_t, Positive_Matrix_With_Prefactor> &parsed_matrices);

void read_mathematica(
  const std::filesystem::path &input_path,
  const std::function<bool(size_t matrix_index)> &should_parse_matrix,
  std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization, size_t &num_matrices,
  std::map<size_t, Positive_Matrix_With_Prefactor> &parsed_matrices)
{
  std::ifstream input_stream(input_path);
  if(!input_stream.good())
    {
      El::RuntimeError("Unable to open input: ", input_path);
    }

  boost::interprocess::file_mapping mapped_file(
    input_path.c_str(), boost::interprocess::read_only);
  boost::interprocess::mapped_region mapped_region(
    mapped_file, boost::interprocess::read_only);

  try
    {
      const char *begin(
        static_cast<const char *>(mapped_region.get_address())),
        *end(begin + mapped_region.get_size());
      parse_SDP(begin, end, should_parse_matrix, objectives, normalization,
                num_matrices, parsed_matrices);
    }
  catch(std::exception &e)
    {
      El::RuntimeError("Error when parsing ", input_path, ": ", e.what());
    }
}
