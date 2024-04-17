#include "Pmp2functions_Parameters.hxx"
#include "pmp_read/pmp_read.hxx"
#include "sdpb_util/assert.hxx"

#include <boost/program_options.hpp>
#include <filesystem>

namespace fs = std::filesystem;
namespace po = boost::program_options;

void write_functions(const fs::path &output_path,
                     const Polynomial_Matrix_Program &pmp);

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      // TODO fix parallel
      if(El::mpi::Size() > 1)
        RUNTIME_ERROR("pmp2functions cannot work in parallel!");

      const Pmp2functions_Parameters parameters(argc, argv);
      if(!parameters.is_valid())
        return 0;

      Environment::set_precision(parameters.precision);

      Timers timers(env, parameters.verbosity);
      const auto pmp
        = read_polynomial_matrix_program(env, parameters.input_file, timers);
      write_functions(parameters.output_path, pmp);
    }
  catch(std::exception &e)
    {
      std::cerr << "Error: " << e.what() << "\n" << std::flush;
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
  catch(...)
    {
      std::cerr << "Unknown Error\n" << std::flush;
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
