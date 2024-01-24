#include "sdpb_util/assert.hxx"

#include <filesystem>
#include <vector>
#include <iostream>
#include <stdexcept>

using namespace std::literals;
namespace fs = std::filesystem;

void parse_command_line(int argc, char **argv, int &precision,
                        std::vector<fs::path> &input_files,
                        fs::path &output_dir)
{
  std::string usage("pvm2sdp [PRECISION] [INPUT]... [OUTPUT_DIR]\n");
  for(int arg = 1; arg < argc; ++arg)
    {
      if((argv[arg] == "-h"s) || argv[arg] == "--help"s)
        {
          std::cerr << usage;
          exit(0);
        }
    }

  if(argc < 4)
    {
      std::cerr << "Wrong number of arguments\n" << usage;
      exit(1);
    }

  std::string precision_string(argv[1]);
  size_t pos;
  try
    {
      precision = (std::stoi(precision_string, &pos));
    }
  catch(std::logic_error &e)
    {
      RUNTIME_ERROR("Invalid precision: '" + precision_string + "'");
    }
  ASSERT(pos == precision_string.size(),
         "Precision has trailing characters: '",
         precision_string.substr(pos) + "'");

  input_files.insert(input_files.end(), argv + 2, argv + argc - 1);
  for(auto &file : input_files)
    {
      ASSERT(fs::exists(file), "Input file does not exist: ", file);
    }
  output_dir = argv[argc - 1];
}
