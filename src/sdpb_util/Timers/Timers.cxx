#include "Timers.hxx"

#include "sdpb_util/Proc_Meminfo.hxx"

namespace
{
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

  void print_meminfo(const std::string &name, bool print_mem_total = false)
  {
    auto prefix = El::BuildString(El::mpi::Rank(), " ", name, " ");

    bool result;
    const auto meminfo = Proc_Meminfo::try_read(result);
    if(!result)
      {
        El::Output(prefix, "cannot parse /proc/meminfo");
        return;
      }

    if(print_mem_total)
      {
        double mem_total_GB
          = static_cast<double>(meminfo.mem_total) / 1024 / 1024 / 1024;
        El::Output(prefix, "MemTotal, GB: ", mem_total_GB);
      }

    double mem_used_GB
      = static_cast<double>(meminfo.mem_used()) / 1024 / 1024 / 1024;
    El::Output(prefix, "MemUsed, GB: ", mem_used_GB);
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
  if(debug && comm_shared_mem.Rank() == 0)
    {
      // Print memory usage for each node (at first rank)
      // MemTotal is constant, thus we print it only once
      // print MemUsed each time
      bool print_mem_total = empty();
      print_meminfo(full_name, print_mem_total);
    }

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
