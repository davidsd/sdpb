#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <system_error>
#include <variant>
#include <vector>

struct Command
{
  using Named_Args_Map = std::map<std::string, std::string>;

  std::string exe;
  std::vector<std::string> args;

  Command(const std::string &exe, const std::vector<std::string> &args,
          const Named_Args_Map &named_args = {});
  Command(const std::string &exe, const Named_Args_Map &named_args);
  [[nodiscard]] static Command split(const std::string &cmd);

  Command &operator+=(const std::string &arg);
  Command &operator+=(const std::vector<std::string> &more_args);
  Command &operator+=(const Named_Args_Map &more_args);
  Command &operator+=(const Command &other);

  template <class T> friend Command operator+(Command lhs, const T &rhs)
  {
    lhs += rhs;
    return lhs;
  }

  std::string to_string() const;
};

std::ostream &operator<<(std::ostream &os, const Command &cmd);

struct Process_Stdio
{
  // NB: separate variant option for nullptr is necessary for constructing bp::process_stdio.
  // Passing FILE* nullptr will lead to segfault.
  using In = std::variant<std::nullptr_t, std::FILE *, std::string>;
  using Out = std::variant<std::nullptr_t, std::FILE *, std::string>;
  using Err = Out;

  In in = stdin;
  Out out = stdout;
  Err err = stderr;
};

int run_command(const Command &command, const Process_Stdio &,
                std::error_code &ec);
int run_command(const Command &command, const Process_Stdio & = {});

std::filesystem::path find_executable(const std::string &path);
