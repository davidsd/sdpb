#include <boost/filesystem.hpp>
#include <vector>

void output_escaped_string(std::ostream &os, const std::string &s)
{
  for(auto &c: s)
    {
      if(c=='"')
        {
          os << "\\\"\\";
        }
      else
        {
          os << c;
        }
    }
}

void write_control(const boost::filesystem::path &output_dir,
                   const int &num_procs,
                   const std::vector<std::string> &command_arguments)
{
  const boost::filesystem::path output_path(
    output_dir / ("control.json"));
  boost::filesystem::ofstream output_stream(output_path);
  output_stream << "{\n  \"num_procs\": " << num_procs
                << ",\n  \"command\": \"";
  for(auto argument(command_arguments.begin()); argument!=command_arguments.end(); ++argument)
    {
      if(argument!=command_arguments.begin())
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
