#include "Pmp2functions_Parameters.hxx"
#include "pmp_read/pmp_read.hxx"
#include "sdpb_util/assert.hxx"

#include <boost/program_options.hpp>
#include <filesystem>

namespace fs = std::filesystem;
namespace po = boost::program_options;

void write_functions(const fs::path &output_path,
                     const std::vector<El::BigFloat> &objectives,
                     const std::vector<El::BigFloat> &normalization,
                     const std::vector<Polynomial_Vector_Matrix> &matrices);

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

      El::gmp::SetPrecision(parameters.precision);
      // El::gmp wants base-2 bits, but boost::multiprecision wants
      // base-10 digits.
      Boost_Float::default_precision(parameters.precision * log(2) / log(10));

      Timers timers(env, parameters.verbosity >= debug);
      const auto pmp
        = read_polynomial_matrix_program(parameters.input_file, timers);
      write_functions(parameters.output_path, pmp.objective, pmp.normalization,
                      pmp.matrices);
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
