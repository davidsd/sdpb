#include "Test_Case_Runner.hxx"
#include "Test_Config.hxx"
#include "process.hxx"
#include "sdpb_util/assert.hxx"

#include <fstream>
#include <boost/algorithm/string/join.hpp>

namespace fs = std::filesystem;

namespace Test_Util
{
  // NB: name should be a valid path relative to test_log_dir
  Test_Case_Runner::Test_Case_Runner(const std::string &name)
      : name(name),
        data_dir(Test_Config::test_data_dir / name),
        output_dir(Test_Config::test_output_dir / name),
        stdout_path(Test_Config::test_log_dir
                    / fs::path(name + ".stdout.log")),
        stderr_path(Test_Config::test_log_dir / fs::path(name + ".stderr.log"))
  {
    fs::remove_all(output_dir);
    fs::remove_all(Test_Config::test_log_dir / name);

    if(stdout_path.has_parent_path())
      fs::create_directories(stdout_path.parent_path());
    ASSERT(fs::is_directory(stdout_path.parent_path()),
           stdout_path.parent_path(),
           " is not a directory! Check file name and permissions.");
  }

  Test_Case_Runner
  Test_Case_Runner::create_nested(const std::string &suffix,
                                  const std::string &separator) const
  {
    return Test_Case_Runner(name + separator + suffix);
  }

  void Test_Case_Runner::run(const Command &command, int required_exit_code,
                             const std::string &required_error_msg) const
  {
    CAPTURE(command);
    CAPTURE(stdout_path);
    CAPTURE(stderr_path);

    // int exit_code = bp::system(command, bp::std_out > stdout_path.string(),
    //                            bp::std_err > stderr_path.string());
    int exit_code = run_command(command, {stdin, stdout_path, stderr_path});
    // NB: We need separate stderr output to process stderr_string.
    // TODO: ideally, we want to redirect bp::std_err to both stdout_path and stderr_path
    // instead of appending stderr to the end ot stdout.
    // Unfortunately, the following simple solution doesn't do what we need,
    // since redirections works only for one path:
    // bp::system(command, (bp::std_out & bp::std_err) > stdout_path, bp::std_err > stderr_path)

    std::ofstream os_stdout(stdout_path, std::ios::app);
    os_stdout << std::endl << "===Command line===" << std::endl;
    os_stdout << command << std::endl;
    os_stdout << "===================" << std::endl;

    const auto stderr_string = [&] {
      std::ifstream is_stderr(stderr_path);
      std::stringstream ss;
      ss << is_stderr.rdbuf();
      return ss.str();
    }();

    os_stdout << std::endl << "===stderr===" << std::endl;
    os_stdout << stderr_string << std::endl;
    os_stdout << "============" << std::endl;

    CAPTURE(exit_code);
    CAPTURE(required_exit_code);
    CAPTURE(stderr_string);
    if(required_exit_code == 0)
      {
        REQUIRE(exit_code == 0);
      }
    else
      {
        REQUIRE(exit_code != 0);
        if(exit_code != required_exit_code)
          {
            // We can get, e.g., code 137 from srun,
            // when MPI process exits with code 1.
            // Test should not fail in this case,
            // so we only print warning.
            WARN("Expected exit code " << exit_code << ", got "
                                       << required_exit_code);
          }
      }
    if(required_exit_code != 0 && !required_error_msg.empty())
      {
        CAPTURE(required_error_msg);
        bool found_error_message
          = stderr_string.find(required_error_msg) != std::string::npos;
        REQUIRE(found_error_message);
      }
  }

  void Test_Case_Runner::mpi_run(const Command &command, const int numProcs,
                                 const int required_exit_code,
                                 const std::string &required_error_msg) const
  {
    auto mpi_command = Test_Config::mpirun;
    mpi_command += {"-n", std::to_string(numProcs)};
    mpi_command += command;
    run(mpi_command, required_exit_code, required_error_msg);
  }

  fs::path Test_Case_Runner::unzip_to_temp_dir(const fs::path &zip_path) const
  {
    const auto temp_dir = output_dir;
    fs::create_directories(temp_dir);
    static int unique_suffix;
    const auto filename
      = zip_path.filename().string() + "." + std::to_string(unique_suffix++);
    auto output_path = temp_dir / filename;
    run({"unzip", {"-o", zip_path, "-d", output_path}});
    return output_path;
  }
}