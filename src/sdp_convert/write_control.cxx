#include "byte_counter.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <vector>

void output_escaped_string(std::ostream &os, const std::string &s)
{
  for(auto &c : s)
    {
      if(c == '"')
        {
          os << "\\\"\\";
        }
      else
        {
          os << c;
        }
    }
}

size_t write_control(const boost::filesystem::path &output_dir,
                     const size_t &num_blocks,
                     const std::vector<std::string> &command_arguments)
{
  const boost::filesystem::path output_path(output_dir / ("control.json"));

  byte_counter counter;
  {
    boost::iostreams::filtering_ostream output_stream;
    output_stream.push(boost::ref(counter));
    // Use gzip with no compression to get a CRC
    output_stream.push(
      boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(0)));
    output_stream.push(boost::iostreams::file_sink(output_path.string()));
    output_stream << "{\n  \"num_blocks\": " << num_blocks
                  << ",\n  \"command\": \"";
    for(auto argument(command_arguments.begin());
        argument != command_arguments.end(); ++argument)
      {
        if(argument != command_arguments.begin())
          {
            output_stream << " ";
          }
        output_escaped_string(output_stream, *argument);
      }
    output_stream << "\"\n}\n";
    if(!output_stream.good())
      {
        throw std::runtime_error("Error when writing to: "
                                 + output_path.string());
      }
  }
  return counter.num_bytes;
}
