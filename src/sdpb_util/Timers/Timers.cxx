#include "Timers.hxx"

#include "sdpb_util/Proc_Meminfo.hxx"

namespace
{
  // Convert bytes to gigabytes
  double to_GB(size_t bytes)
  {
    return static_cast<double>(bytes) / 1024 / 1024 / 1024;
  }

  // TODO print_statm() is currently unused

  // /proc/self/statm displays the following quantities:
  // size resident shared text lib data dt
  // Each quantity is measured in pages;
  // usually page size=4096B (you can check it by calling "getconf PAGESIZE")
  void print_statm(const std::string &prefix)
  {
    std::ifstream stat_file("/proc/self/statm");
    if(stat_file.good())
      {
        std::string stats;
        std::getline(stat_file, stats);
        El::Output(prefix, stats);
      }
  }
}

Timers::Timers(bool debug) : debug(debug)
{
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &comm_shared_mem.comm);
}
Timer &Timers::add_and_start(const std::string &name)
{
  std::string full_name = prefix + name;

  if(debug)
    print_meminfo(full_name);

  emplace_back(full_name, Timer());
  return back().second;
}
void Timers::write_profile(const std::filesystem::path &path) const
{
  if(!path.parent_path().empty())
    std::filesystem::create_directories(path.parent_path());
  std::ofstream f(path);

  f << "{" << '\n';
  for(auto it(begin()); it != end();)
    {
      f << "    {\"" << it->first << "\", " << it->second << "}";
      ++it;
      if(it != end())
        {
          f << ",";
        }
      f << '\n';
    }
  f << "}" << '\n';

  if(!f.good())
    {
      throw std::runtime_error("Error when writing to: " + path.string());
    }
}
int64_t Timers::elapsed_milliseconds(const std::string &s) const
{
  auto iter(std::find_if(rbegin(), rend(),
                         [&s](const std::pair<std::string, Timer> &timer) {
                           return timer.first == s;
                         }));
  if(iter == rend())
    {
      throw std::runtime_error("Could not find timing for " + s);
    }
  return iter->second.elapsed_milliseconds();
}

void Timers::print_max_mem_used() const
{
  if(max_mem_used > 0 && !max_mem_used_name.empty())
    {
      El::Output(El::mpi::Rank(), " max MemUsed: ", to_GB(max_mem_used),
                 " GB at \"", max_mem_used_name, "\"");
    }
}

void Timers::print_meminfo(const std::string &name)
{
  // Print data from /proc/meminfo only for a first rank of each node
  if(comm_shared_mem.Rank() != 0)
    return;

  auto prefix = El::BuildString(El::mpi::Rank(), " ", name, " ");

  // Print memory usage for the current node (from the first rank).
  // If we cannot parse /proc/meminfo, then simply print timer name.

  if(!can_read_meminfo)
    {
      El::Output(prefix);
      return;
    }

  bool result;
  constexpr bool print_error_msg = true;
  const auto meminfo = Proc_Meminfo::try_read(result, print_error_msg);
  if(!result)
    {
      can_read_meminfo = false;
      El::Output("Printing RAM usage will be disabled.");
      El::Output(prefix);
      return;
    }

  // MemTotal is constant, thus we print it only once, when adding first timer
  if(empty())
    {
      El::Output(prefix, "--- MemTotal: ", to_GB(meminfo.mem_total), " GB");
    }

  //Print MemUsed each time
  El::Output(prefix, "--- MemUsed: ", to_GB(meminfo.mem_used()), " GB");

  // Update max MemUsed info
  if(meminfo.mem_used() > max_mem_used)
    {
      max_mem_used = meminfo.mem_used();
      max_mem_used_name = name;
    }
}
