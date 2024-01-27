#include <filesystem>
#include <vector>
#include <iostream>
#include <stdexcept>
#include "pmp2sdp/Block_File_Format.hxx"

using namespace std::literals;
namespace fs = std::filesystem;

void parse_command_line(int argc, char **argv,
                        Block_File_Format &output_format, int &precision,
                        std::vector<fs::path> &input_files,
                        fs::path &output_dir,
                        std::vector<std::string> &command_arguments)
{
  std::string usage("pvm2sdp [FORMAT] PRECISION INPUT... OUTPUT\n"
                    "FORMAT (optional): output format, bin (default) or json\n"
                    "PRECISION: binary precision\n"
                    "INPUT: input .xml or .nsv files\n"
                    "OUTPUT: output sdp directory");
  for(int arg = 1; arg < argc; ++arg)
    {
      if((argv[arg] == "-h"s) || argv[arg] == "--help"s)
        {
          std::cerr << usage;
          exit(0);
        }
    }

  for(int arg(0); arg != argc; ++arg)
    {
      command_arguments.emplace_back(argv[arg]);
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
      RUNTIME_ERROR("Invalid precision: '" + precision_string + "'");
    }
  if(pos != precision_string.size())
    {
      RUNTIME_ERROR("Precision has trailing characters: '"
                    + precision_string.substr(pos) + "'");
    }
  curr_arg_pos++;

  input_files.insert(input_files.end(), argv + curr_arg_pos, argv + argc - 1);
  for(auto &file : input_files)
    {
      ASSERT(fs::exists(file), "Input file does not exist: ", file);
    }
  output_dir = argv[argc - 1];
}
