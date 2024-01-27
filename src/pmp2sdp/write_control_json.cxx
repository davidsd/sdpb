#include <iostream>
#include <string>
#include <vector>
#include "sdpb_util/ostream/set_stream_precision.hxx"

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

void write_control_json(std::ostream &output_stream, const size_t &num_blocks,
                        const std::vector<std::string> &command_arguments)

{
  set_stream_precision(output_stream);
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
}
