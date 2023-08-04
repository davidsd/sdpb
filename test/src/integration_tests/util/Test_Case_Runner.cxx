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

  int Test_Case_Runner::run(const std::string &command) const
  {
    namespace bp = boost::process;

    boost::filesystem::ofstream os_stdout(stdout_path, std::ios::app);
    boost::filesystem::ofstream os_stderr(stderr_path, std::ios::app);

    // write command before stdout
    os_stdout << command << std::endl;

    bp::ipstream stdout_pipe;
    bp::ipstream stderr_pipe;
    int result = bp::system(command, bp::std_out > stdout_pipe,
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

    return result;
  }

  int Test_Case_Runner::mpi_run(const std::string &command, int numProcs) const
  {
    auto mpi_command
      = build_command_line(build_mpirun_prefix(numProcs), command);
    return run(mpi_command);
  }

  int Test_Case_Runner::run(const std::vector<std::string> &args,
                            const Named_Args_Map &named_args) const
  {
    auto args_string = boost::algorithm::join(args, " ");
    auto named_args_string = build_string_from_named_args(named_args);
    auto command = build_command_line(args_string, named_args_string);
    return run(command);
  }

  int Test_Case_Runner::mpi_run(const std::vector<std::string> &args,
                                const Named_Args_Map &named_args,
                                int numProcs) const
  {
    std::vector<std::string> args_with_mpi(args);
    args_with_mpi.insert(args_with_mpi.begin(), build_mpirun_prefix(numProcs));
    return run(args_with_mpi, named_args);
  }

  // Note that calling diff on big files is slow and generates long output.
  // Thus, instead of diff, we call cmp - it exits at first difference.
  int Test_Case_Runner::diff(const boost::filesystem::path &x,
                             const boost::filesystem::path &y) const
  {
    if(is_directory(x) != is_directory(y))
      {
        WARN("Cannot compare file and directory: " << x << " " << y << "");
        return 1;
      }

    if(is_directory(x))
      {
        if(!is_directory(y))
          {
            WARN(y << " should be directory");
            return 1;
          }

        std::vector<boost::filesystem::path> x_children;
        std::vector<boost::filesystem::path> y_children;

        for(const auto &p : boost::make_iterator_range(
              boost::filesystem::directory_iterator(x), {}))
          {
            x_children.push_back(p);
          }
        for(const auto &p : boost::make_iterator_range(
              boost::filesystem::directory_iterator(y), {}))
          {
            y_children.push_back(p);
          }
        if(x_children.size() != y_children.size())
          {
            WARN(x << " has " << x_children.size() << " entries, " << y
                   << " has " << y_children.size() << " entries.");
            return 1;
          }
        for(const auto &x_child : x_children)
          {
            auto y_child = y / boost::filesystem::relative(x_child, x);
            int res = diff(x_child, y_child);
            if(res != 0)
              return res;
          }
        return 0;
      }

    // diff can be slow for large files, thus we call cmp
    auto command = build_command_line("cmp --print-bytes", x, y);
    return run(command);
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
    CAPTURE(unzip_command);
    REQUIRE(run(unzip_command) == 0);
    return output_path;
  }

  bool Test_Case_Runner::stderr_contains_substring(
    const std::string &substring) const
  {
    boost::filesystem::ifstream is(stderr_path);

    std::string line;
    while(std::getline(is, line))
      {
        if(line.find(substring) != std::string::npos)
          return true;
      }
    return false;
  }
}