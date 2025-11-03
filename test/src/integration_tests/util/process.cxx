#include "process.hxx"

#include "sdpb_util/assert.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"

#include <type_traits>
#include <boost/version.hpp>
#include <boost/program_options/parsers.hpp>

// There are two versions for Boost.Process, V1 and V2.
// V2 was introduced in Boost 1.80 and became default in Boot 1.89.
// We use V2 when possible.
#if BOOST_VERSION < 108000
#define USE_BOOST_PROCESS_V1
#include <boost/process.hpp>
namespace bp = boost::process;
#elif BOOST_VERSION < 108900
#include <boost/process/v2.hpp>
namespace bp = boost::process::v2;
#else
#include <boost/process.hpp>
namespace bp = boost::process;
#endif

// struct Command

Command::Command(const std::string &exe, const std::vector<std::string> &args,
                 const Named_Args_Map &named_args)
    : exe(exe), args(args)
{
  *this += named_args;
}
Command::Command(const std::string &exe, const Named_Args_Map &named_args)
    : Command(exe, {}, named_args)
{}

Command Command::split(const std::string &cmd)
{
  const auto cmd_with_args = boost::program_options::split_unix(cmd);
  ASSERT(!cmd_with_args.empty(), "Cannot split cmd=", cmd);
  auto begin = cmd_with_args.begin();
  const auto end = cmd_with_args.end();
  const auto exe = *begin++;
  const std::vector args(begin, end);
  return Command(exe, args);
}
Command &Command::operator+=(const std::string &arg)
{
  if(!arg.empty())
    args.push_back(arg);
  return *this;
}
Command &Command::operator+=(const std::vector<std::string> &more_args)
{
  args.reserve(args.size() + more_args.size());
  for(auto &arg : more_args)
    *this += arg;
  return *this;
}
Command &Command::operator+=(const Named_Args_Map &more_args)
{
  args.reserve(args.size() + more_args.size() * 2);
  for(const auto &[key, value] : more_args)
    {
      *this += key;
      *this += value;
    }
  return *this;
}
Command &Command::operator+=(const Command &other)
{
  *this += other.exe;
  return *this += other.args;
}
std::string Command::to_string() const
{
  std::ostringstream os;
  os << *this;
  return os.str();
}
std::ostream &operator<<(std::ostream &os, const Command &cmd)
{
  os << cmd.exe;
  for(const auto &arg : cmd.args)
    {
      os << ' ' << arg;
    }
  return os;
}

// struct Command

template <typename T, typename Variant> struct is_in_variant;

template <typename T, typename... Ts>
struct is_in_variant<T, std::variant<Ts...>>
    : std::disjunction<std::is_same<T, Ts>...>
{};

template <typename T, typename Variant>
inline constexpr bool is_in_variant_v = is_in_variant<T, Variant>::value;

// Helper structure for Process_Stdio.
// TODO: shall we use it in API instead of Process_Stdio (which uses std::variant)?
// This would allow us to replace std::visit with template functions,
// but then we'd have to move implementation to header file
template <class TIn, class TOut, class TErr,
          typename
          = std::enable_if_t<is_in_variant_v<TIn, Process_Stdio::In>
                               && is_in_variant_v<TOut, Process_Stdio::Out>
                               && is_in_variant_v<TErr, Process_Stdio::Err>,
                             int>>
struct TProcess_Stdio
{
  TIn in;
  TOut out;
  TErr err;

  TProcess_Stdio(const TIn &in, const TOut &out, const TErr &err)
      : in(in), out(out), err(err)
  {}
};

#ifdef USE_BOOST_PROCESS_V1
template <class T> auto to_bp_redirect_arg(const T &value)
{
  return value;
}
template <> auto to_bp_redirect_arg<std::nullptr_t>(const std::nullptr_t &)
{
  return bp::null;
}
#endif

template <class TIn, class TOut, class TErr>
int run_command(const Command &command,
                const TProcess_Stdio<TIn, TOut, TErr> &proc_stdio,
                std::error_code &ec)
{
  try
    {
#ifdef USE_BOOST_PROCESS_V1
      const auto in = bp::std_in < to_bp_redirect_arg(proc_stdio.in);
      const auto out = bp::std_out > to_bp_redirect_arg(proc_stdio.out);
      const auto err = bp::std_err > to_bp_redirect_arg(proc_stdio.err);
      return bp::system(command.to_string(), in, out, err, ec);
#else
      boost::asio::io_context ctx;
      bp::process_stdio io{bp::detail::process_input_binding(proc_stdio.in),
                           bp::detail::process_output_binding(proc_stdio.out),
                           bp::detail::process_error_binding(proc_stdio.err)};
      boost::system::error_code local_ec;
      const auto exe = find_executable(command.exe).string();
      ASSERT(!exe.empty(), "Cannot find executable: ", command.exe);
      const auto &args = command.args;
      return bp::execute(bp::process(ctx, exe, args, io), local_ec);
      ec.assign(local_ec.value(), local_ec.category());
#endif
    }
  catch(std::exception &ex)
    {
      std::ostringstream os;
      os << "Failed to run command: exe=" << command.exe
         << ", args=" << command.args << ":\n"
         << ex.what();
      RUNTIME_ERROR(os.str());
    }
}

int run_command(const Command &command, const Process_Stdio &proc_stdio,
                std::error_code &ec)
{
  return std::visit(
    [&](auto &&in) {
      return std::visit(
        [&](auto &&out) {
          return std::visit(
            [&](auto &&err) {
              return run_command(command, TProcess_Stdio(in, out, err), ec);
            },
            proc_stdio.err);
        },
        proc_stdio.out);
    },
    proc_stdio.in);
}

int run_command(const Command &command, const Process_Stdio &proc_stdio)
{
  std::error_code ec;
  const int exit_code = run_command(command, proc_stdio, ec);
  ASSERT(!ec, ec.message());
  return exit_code;
}

std::filesystem::path find_executable(const std::string &filename)
{
  // Boost.Process returns empty string for paths with directory separators.
  // We still want an executable, so we simply return the path.
  // TODO: check also that the file exists and is executable?
  if(std::filesystem::path(filename).has_parent_path())
    return filename;
#ifdef USE_BOOST_PROCESS_V1
  return bp::search_path(filename).string();
#else
  return bp::environment::find_executable(filename).string();
#endif
}
