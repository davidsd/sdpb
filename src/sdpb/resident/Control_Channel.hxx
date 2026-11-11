#pragma once

#include "sdpb_util/Environment.hxx"

#include <filesystem>
#include <optional>
#include <string>

// File mailbox between a driver and a resident SDPB (sdpb --controlDir).
//
// The driver writes the commands <dir>/cmd.<seq>.json for seq = 1, 2, ...
// Each file must appear at once: write it under another name in <dir>, then
// rename it. SDPB answers each command with <dir>/res.<seq>.json (also
// renamed into place). SDPB writes <dir>/ready.json once the SDP is loaded
// and <dir>/exit.json when it stops, and touches <dir>/alive every ~10 s
// while it waits. A file <dir>/stop makes it stop at the next wait.
//
// Only rank 0 touches the files. The other ranks wait in MPI_Ibarrier and
// sleep, so that waiting ranks leave the CPUs to the driver.
class Control_Channel
{
public:
  enum class Stop_Reason
  {
    none,
    idle_timeout,
    sigterm,
    stop_file
  };

  Control_Channel(const std::filesystem::path &dir,
                  int64_t idle_timeout_seconds);

  // Wait for the next command. Collective. Returns the text of the command on
  // all ranks, or nullopt on idle timeout, SIGTERM or a stop file (see
  // stop_reason()).
  std::optional<std::string> next_command(const Environment &env);

  [[nodiscard]] size_t seq() const { return current_seq; }
  [[nodiscard]] Stop_Reason stop_reason() const { return reason; }
  [[nodiscard]] const std::filesystem::path &directory() const { return dir; }

  // Rank 0 only: write res.<seq>.json for the current command
  void write_result(const std::string &json) const;
  // Rank 0 only: write <dir>/<name> atomically
  void write_file(const std::string &name, const std::string &text) const;

private:
  std::filesystem::path dir;
  int64_t idle_timeout_seconds;
  // Parent process at start; a waiting rank whose parent changed exits
  int parent_pid;
  [[nodiscard]] bool orphaned() const;
  size_t current_seq = 0;
  Stop_Reason reason = Stop_Reason::none;
};

std::string to_string(Control_Channel::Stop_Reason reason);
