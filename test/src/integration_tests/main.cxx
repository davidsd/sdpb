#include "common.hxx"
#include "sdpb_util/assert.hxx"

#include <filesystem>
#include <boost/process.hpp>
#include <iostream>

#ifndef CATCH_AMALGAMATED_CUSTOM_MAIN
#error "To override main, pass '-D CATCH_AMALGAMATED_CUSTOM_MAIN' to compiler"
#endif

namespace bp = boost::process;
namespace fs = std::filesystem;

namespace
{
  // Change working directory to sdpb root
  // Allow to run the executable from sdpb/ or sdpb/build/ directory,
  // otherwise fail
  void set_working_directory()
  {
    // if we find this file, we know where we are
    const auto sdp_zip = "test/data/pmp2sdp/xml/sdp_orig.zip";
    for(const fs::path prefix : {".", ".."})
      {
        if(exists(prefix / sdp_zip))
          {
            fs::current_path(prefix);
            return;
          }
      }

    RUNTIME_ERROR("Cannot find '", sdp_zip,
                  "'. Please run test executable from root sdpb/ directory "
                  "or sdpb/build/ directory");
  }
}

// main() with custom command line option --mpirun,
// based on
// https://github.com/catchorg/Catch2/blob/devel/docs/own-main.md#adding-your-own-command-line-options
int main(int argc, char *argv[])
{
  Catch::Session session; // There must be exactly one instance

  set_working_directory();

  // Build a new command line parser on top of Catch2's
  using namespace Catch::Clara;
  std::string mpirun = "mpirun";
  // bind mpirun_command variable to a new option --mpirun
  auto cli = session.cli()
             | Opt(mpirun, "mpirun")["--mpirun"](
               "mpirun command, e.g. --mpirun=srun or "
               "--mpirun=\"mpirun --mca btl vader,openib,self\"."
               "By default, --mpirun=mpirun");

  // Now pass the new composite back to Catch2 so it uses that
  session.cli(cli);

  // Let Catch2 (using Clara) parse the command line
  int returnCode = session.applyCommandLine(argc, argv);
  if(returnCode != 0) // Indicates a command line error
    return returnCode;

  // TODO process arguments like --mpirun="mpirun --mca btl vader,openib,self"
  Test_Config::mpirun = mpirun;

  // Check if we can run mpirun
  auto mpirun_help_command = mpirun + " --help";
  std::error_code ec;
  int check_mpi = bp::system(mpirun_help_command, bp::std_out > bp::null, ec);
  if(check_mpi != 0)
    {
      std::cout << "Failed to run MPI: " << mpirun_help_command << std::endl;
      if(ec.value() != 0)
        std::cout << "Error code: " << ec << std::endl;
      else
        std::cout << "Exit code: " << check_mpi << std::endl;
      std::cout << "Please specify mpirun command, e.g.:" << std::endl
                << "./build/integration_tests --mpirun=srun" << std::endl
                << "./build/integration_tests"
                   " --mpirun=\"mpirun --mca btl vader,openib,self\""
                << std::endl;
      return 1;
    }

  // clear output and log directory
  fs::remove_all(Test_Config::test_output_dir);
  fs::remove_all(Test_Config::test_log_dir);

  return session.run();
}
