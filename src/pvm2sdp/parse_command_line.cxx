#include <boost/filesystem.hpp>
#include <vector>
#include <iostream>
#include <stdexcept>
#include "../sdp_convert/Block_File_Format.hxx"

using namespace std::literals;

void parse_command_line(int argc, char **argv,
                        Block_File_Format &output_format, int &precision,
                        std::vector<boost::filesystem::path> &input_files,
                        boost::filesystem::path &output_dir)
{
  std::string usage("pvm2sdp [FORMAT] PRECISION INPUT... OUTPUT\n"
                    "FORMAT (optional): output format, bin (default) or json\n"
                    "PRECISION: binary precision\n"
                    "INPUT: input .xml or .nsv files\n"
                    "OUTPUT: output sdp.zip file");
  for(int arg = 1; arg < argc; ++arg)
    {
      if((argv[arg] == "-h"s) || argv[arg] == "--help"s)
        {
          std::cerr << usage;
          exit(0);
        }
    }

  size_t curr_arg_pos = 1;
  std::string first_arg = argv[curr_arg_pos];
  if(first_arg == "bin" || first_arg == "json")
    {
      output_format = first_arg == "bin" ? bin : json;
      curr_arg_pos = 2;
    }

  if(argc - curr_arg_pos < 3)
    {
      std::cerr << "Wrong number of arguments\n" << usage;
      exit(1);
    }

  std::string precision_string(argv[curr_arg_pos]);
  size_t pos;
  try
    {
      precision = (std::stoi(precision_string, &pos));
    }
  catch(std::logic_error &e)
    {
      throw std::runtime_error("Invalid precision: '" + precision_string
                               + "'");
    }
  if(pos != precision_string.size())
    {
      throw std::runtime_error("Precision has trailing characters: '"
                               + precision_string.substr(pos) + "'");
    }
  curr_arg_pos++;

  input_files.insert(input_files.end(), argv + curr_arg_pos, argv + argc - 1);
  for(auto &file : input_files)
    {
      if(!boost::filesystem::exists(file))
        {
          throw std::runtime_error("Input file '" + file.string()
                                   + "' does not exist");
        }
    }
  output_dir = argv[argc - 1];
}
