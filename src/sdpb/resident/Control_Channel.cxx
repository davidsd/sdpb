#include "Control_Channel.hxx"

#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <unistd.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <thread>

namespace fs = std::filesystem;

namespace
{
  std::string read_whole_file(const fs::path &path)
  {
    std::ifstream in(path, std::ios::binary);
    ASSERT(in.good(), "Cannot open ", path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  void broadcast_string(std::string &s)
  {
    size_t size = s.size();
    El::mpi::Broadcast(size, 0, El::mpi::COMM_WORLD);
    s.resize(size);
    if(size > 0)
      MPI_Bcast(s.data(), size, MPI_CHAR, 0, El::mpi::COMM_WORLD.comm);
  }
}

std::string to_string(const Control_Channel::Stop_Reason reason)
{
  switch(reason)
    {
    case Control_Channel::Stop_Reason::none: return "none";
    case Control_Channel::Stop_Reason::idle_timeout: return "idle timeout";
    case Control_Channel::Stop_Reason::sigterm: return "SIGTERM";
    case Control_Channel::Stop_Reason::stop_file: return "stop file";
    }
  return "unknown";
}

Control_Channel::Control_Channel(const fs::path &dir,
                                 const int64_t idle_timeout_seconds)
    : dir(dir), idle_timeout_seconds(idle_timeout_seconds),
      parent_pid(getppid())
{
  if(El::mpi::Rank() == 0)
    {
      fs::create_directories(dir);
      // Results and status files of an earlier process would be mistaken for
      // ours; commands written in advance are kept
      for(const auto &entry : fs::directory_iterator(dir))
        {
          const auto name = entry.path().filename().string();
          if(name.rfind("res.", 0) == 0 || name == "ready.json"
             || name == "exit.json" || name == "stop" || name == "alive")
            fs::remove(entry.path());
        }
    }
  El::mpi::Barrier(El::mpi::COMM_WORLD);
}

bool Control_Channel::orphaned() const
{
  // The launcher (mpirun's daemon, slurmstepd) is gone: nobody will send
  // commands or clean up this process
  return getppid() != parent_pid;
}

std::optional<std::string> Control_Channel::next_command(const Environment &env)
{
  using Clock = std::chrono::steady_clock;
  const size_t seq = current_seq + 1;
  std::string text;
  int status = 0;
  if(El::mpi::Rank() == 0)
    {
      const fs::path path = dir / ("cmd." + std::to_string(seq) + ".json");
      const auto start = Clock::now();
      auto last_alive = start - std::chrono::hours(1);
      std::chrono::microseconds sleep(100);
      while(true)
        {
          if(fs::exists(path))
            {
              text = read_whole_file(path);
              break;
            }
          if(env.sigterm_received())
            {
              status = static_cast<int>(Stop_Reason::sigterm);
              break;
            }
          if(orphaned())
            {
              // The other ranks may be gone too: no MPI from here
              write_file("exit.json", "{\"status\": \"exit\", \"reason\": "
                                      "\"launcher exited\", \"pid\": "
                                        + std::to_string(getpid()) + "}");
              std::_Exit(1);
            }
          if(fs::exists(dir / "stop"))
            {
              status = static_cast<int>(Stop_Reason::stop_file);
              break;
            }
          const auto now = Clock::now();
          if(idle_timeout_seconds > 0
             && now - start >= std::chrono::seconds(idle_timeout_seconds))
            {
              status = static_cast<int>(Stop_Reason::idle_timeout);
              break;
            }
          if(now - last_alive >= std::chrono::seconds(10))
            {
              std::ofstream(dir / "alive") << seq << "\n";
              last_alive = now;
            }
          std::this_thread::sleep_for(sleep);
          // Poll fast at first (a driver answering at once), then every 20 ms
          sleep = std::min(sleep * 2, std::chrono::microseconds(20000));
        }
    }

  // Wait for rank 0 without spinning
  MPI_Request request;
  MPI_Ibarrier(El::mpi::COMM_WORLD.comm, &request);
  if(El::mpi::Rank() == 0)
    {
      MPI_Wait(&request, MPI_STATUS_IGNORE);
    }
  else
    {
      int done = 0;
      while(!done)
        {
          MPI_Test(&request, &done, MPI_STATUS_IGNORE);
          if(!done)
            {
              if(orphaned())
                std::_Exit(1);
              std::this_thread::sleep_for(std::chrono::microseconds(500));
            }
        }
    }

  El::mpi::Broadcast(status, 0, El::mpi::COMM_WORLD);
  if(status != 0)
    {
      reason = static_cast<Stop_Reason>(status);
      return std::nullopt;
    }
  broadcast_string(text);
  current_seq = seq;
  return text;
}

void Control_Channel::write_file(const std::string &name,
                                 const std::string &text) const
{
  ASSERT_EQUAL(El::mpi::Rank(), 0);
  const fs::path tmp = dir / (".tmp." + name);
  {
    std::ofstream out(tmp, std::ios::binary);
    out << text;
    ASSERT(out.good(), "Cannot write ", tmp);
  }
  fs::rename(tmp, dir / name);
}

void Control_Channel::write_result(const std::string &json) const
{
  write_file("res." + std::to_string(current_seq) + ".json", json);
}
