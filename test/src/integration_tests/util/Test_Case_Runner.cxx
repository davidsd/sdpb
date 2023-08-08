#include "Test_Case_Runner.hxx"
#include "Test_Config.hxx"

namespace
{
  // concatenate args with " " separator
  inline void build_args_stream(std::ostringstream &) {}
  template <typename T, typename... ArgPack>
  void build_args_stream(std::ostringstream &os, const T &item,
                         const ArgPack &...args)
  {
    os << item << " ";
    build_args_stream(os, args...);
  }

  // concatenate args with " " separator
  inline std::string build_args_string(const std::string &arg) { return arg; }
  template <typename... ArgPack>
  std::string build_args_string(const ArgPack &...args)
  {
    std::ostringstream os;
    build_args_stream(os, args...);
    return os.str();
  }

  template <typename... Ts> std::string build_command_line(const Ts &...args)
  {
    return build_args_string(args...);
  }

  // mpirun -n 2
  std::string build_mpirun_prefix(int numProcs)
  {
    return build_command_line(Test_Config::mpirun, "-n",
                              std::to_string(numProcs));
  }

  std::string build_string_from_named_args(
    const Test_Util::Test_Case_Runner::Named_Args_Map &named_args)
  {
    std::stringstream ss;
    for(const auto &[key, value] : named_args)
      {
        ss << " " << key;
        if(!value.empty())
          ss << "=" << value;
      }
    return ss.str();
  }
}

namespace Test_Util
{
  // NB: name should be a valid path relative to test_log_dir
  Test_Case_Runner::Test_Case_Runner(const std::string &name)
      : name(name), data_dir(Test_Config::test_data_dir / name),
        output_dir(Test_Config::test_output_dir / name),
        stdout_path(Test_Config::test_log_dir
                    / boost::filesystem::path(name + ".stdout.log")),
        stderr_path(Test_Config::test_log_dir
                    / boost::filesystem::path(name + ".stderr.log"))
  {
    boost::filesystem::remove_all(output_dir);
    boost::filesystem::remove_all(Test_Config::test_log_dir / name);

    boost::filesystem::create_directories(stdout_path.parent_path());
    if(!boost::filesystem::is_directory(stdout_path.parent_path()))
      {
        throw std::runtime_error(
          stdout_path.parent_path().string()
          + " is not a directory! Check file name and permissions.");
      }
  }

  Test_Case_Runner
  Test_Case_Runner::create_nested(const std::string &suffix,
                                  const std::string &separator) const
  {
    return Test_Case_Runner(name + separator + suffix);
  }

  void
  Test_Case_Runner::run(const std::string &command, int required_exit_code,
                        const std::string &required_error_msg) const
  {
    namespace bp = boost::process;

    CAPTURE(command);
    CAPTURE(stdout_path);
    CAPTURE(stderr_path);

    boost::filesystem::ofstream os_stdout(stdout_path, std::ios::app);
    boost::filesystem::ofstream os_stderr(stderr_path, std::ios::app);

    // write command before stdout
    os_stdout << command << std::endl;

    bp::ipstream stdout_pipe;
    bp::ipstream stderr_pipe;
    int exit_code = bp::system(command, bp::std_out > stdout_pipe,
                               bp::std_err > stderr_pipe);

    // write stdout to file
    os_stdout << stdout_pipe.rdbuf();

    // write stderr to both stderr and stdout files
    // we cannot call rdbuf() twice, thus we copy it to stderr_string
    std::stringstream ss;
    ss << stderr_pipe.rdbuf();
    auto stderr_string = ss.str();

    os_stdout << stderr_string;
    os_stderr << stderr_string;

    REQUIRE(exit_code == required_exit_code);
    if(required_exit_code != 0 && !required_error_msg.empty())
      {
        CAPTURE(required_error_msg);
        bool found_error_message
          = stderr_string.find(required_error_msg) != std::string::npos;
        REQUIRE(found_error_message);
      }
  }

  void Test_Case_Runner::mpi_run(const std::string &command, int numProcs,
                                 int required_exit_code,
                                 const std::string &required_error_msg) const
  {
    auto mpi_command
      = build_command_line(build_mpirun_prefix(numProcs), command);
    run(mpi_command, required_exit_code, required_error_msg);
  }

  void Test_Case_Runner::run(const std::vector<std::string> &args,
                             const Named_Args_Map &named_args,
                             int required_exit_code,
                             const std::string &required_error_msg) const
  {
    auto args_string = boost::algorithm::join(args, " ");
    auto named_args_string = build_string_from_named_args(named_args);
    auto command = build_command_line(args_string, named_args_string);
    run(command, required_exit_code, required_error_msg);
  }

  void Test_Case_Runner::mpi_run(const std::vector<std::string> &args,
                                 const Named_Args_Map &named_args,
                                 int numProcs, int required_exit_code,
                                 const std::string &required_error_msg) const
  {
    std::vector<std::string> args_with_mpi(args);
    args_with_mpi.insert(args_with_mpi.begin(), build_mpirun_prefix(numProcs));
    run(args_with_mpi, named_args, required_exit_code,
        required_error_msg);
  }

  boost::filesystem::path Test_Case_Runner::unzip_to_temp_dir(
    const boost::filesystem::path &zip_path) const
  {
    auto temp_dir = output_dir;
    boost::filesystem::create_directories(temp_dir);
    auto filename = zip_path.filename().string() + "."
                    + boost::filesystem::unique_path().string();
    auto output_path = temp_dir / filename;
    auto unzip = boost::process::search_path("unzip");
    if(unzip.empty())
      FAIL("Cannot find unzip");

    // control.json may differ by "command" field
    // thus we exclude this file from comparison
    auto unzip_command
      = build_command_line("unzip -o", zip_path, "-d", output_path);
    run(unzip_command);
    return output_path;
  }
}