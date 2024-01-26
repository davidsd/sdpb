#include "Output_SDP/Output_SDP.hxx"
#include "Pmp2sdp_Parameters/Pmp2sdp_Parameters.hxx"
#include "Dual_Constraint_Group.hxx"
#include "write_sdp.hxx"
#include "pmp_read/pmp_read.hxx"
#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <string>
#include <boost/program_options.hpp>
#include <filesystem>

namespace fs = std::filesystem;
namespace po = boost::program_options;

void convert(const Polynomial_Matrix_Program &pmp,
             El::BigFloat &objective_const,
             std::vector<El::BigFloat> &dual_objective_b,
             std::vector<Dual_Constraint_Group> &dual_constraint_groups,
             Timers &timers);

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      Pmp2sdp_Parameters parameters(argc, argv);
      if(!parameters.is_valid())
        {
          return 0;
        }

      El::gmp::SetPrecision(parameters.precision);
      // El::gmp wants base-2 bits, but boost::multiprecision wants
      // base-10 digits.
      Boost_Float::default_precision(parameters.precision * log(2) / log(10));

      const bool debug = parameters.verbosity >= Verbosity::debug;
      Timers timers(env, debug);
      Scoped_Timer timer(timers, "pmp2sdp");

      auto pmp
        = read_polynomial_matrix_program(env, parameters.input_file, timers);

      Output_SDP sdp(pmp, parameters.command_arguments, timers);
      write_sdp(parameters.output_path, sdp, parameters.output_format, timers,
                debug);
      if(El::mpi::Rank() == 0)
        {
          El::Output("Processed ", sdp.num_blocks, " SDP blocks in ",
                     (double)timer.timer().elapsed_milliseconds() / 1000,
                     " seconds, output: ", parameters.output_path.string());
        }

      if(debug)
        {
          timers.write_profile(parameters.output_path.string()
                               + ".profiling/profiling."
                               + std::to_string(El::mpi::Rank()));
        }
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
