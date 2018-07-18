#include <boost/filesystem.hpp>
#include <vector>
#include <iostream>

using namespace std::literals;

int main(int argc, char **argv)
{
  int result(0);
  try
    {
      std::string usage("pvm2sdp [PRECISION] [INPUT]... [OUTPUT_DIR]\n");
      if(argc < 4)
        {
          std::cerr << "Wrong number of arguments\n" << usage;
          return 1;
        }

      for(int arg = 1; arg < argc; ++arg)
        {
          if(argv[arg] == "-h"s || argv[arg] == "--help"s)
            {
              std::cerr << usage;
              return 1;
            }
        }

      int precision([&]() {
        std::string precision_string(argv[1]);
        size_t pos;
        int precision;
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
        return precision;
      }());
      std::vector<boost::filesystem::path> input_files(argv + 2,
                                                       argv + argc - 2);
      boost::filesystem::path output_dir(argv[argc - 1]);
      if(boost::filesystem::exists(output_dir)
         && !boost::filesystem::is_directory(output_dir))
        {
          throw std::runtime_error("Output directory '" + output_dir.string()
                                   + "' exists but is not a directory.");
        }
    }
  catch(std::runtime_error &e)
    {
      std::cerr << "Error: " << e.what() << "\n" << std::flush;
      result = 1;
    }
  catch(...)
    {
      std::cerr << "Unknown Error\n" << std::flush;
      result = 1;
    }
  return result;
}
