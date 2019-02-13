#include "../Input_Parser.hxx"

#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>

std::vector<char>::const_iterator
parse_SDP(const std::vector<char>::const_iterator &begin,
          const std::vector<char>::const_iterator &end,
          std::vector<El::BigFloat> &objectives,
          std::vector<El::BigFloat> &normalization,
          std::vector<Positive_Matrix_With_Prefactor> &matrices);

void read_mathematica(const boost::filesystem::path &input_path,
                      std::vector<El::BigFloat> &objectives,
                      std::vector<El::BigFloat> &normalization,
                      std::vector<Positive_Matrix_With_Prefactor> &matrices)

{
  boost::filesystem::ifstream input_stream(input_path);
  if(!input_stream.good())
    {
      throw std::runtime_error("Unable to open input: " + input_path.string());
    }
  std::vector<char> input_string;
  input_string.resize(boost::filesystem::file_size(input_path));
  input_stream.read(input_string.data(), input_string.size());

  parse_SDP(input_string.cbegin(), input_string.cend(), objectives,
            normalization, matrices);
}
