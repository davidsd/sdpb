#include "Timers.hxx"

#include "sdpb_util/Proc_Meminfo.hxx"
#include "sdpb_util/assert.hxx"

namespace
{
  // Convert bytes to gigabytes
  double to_GB(const size_t bytes)
  {
    return static_cast<double>(bytes) / 1024 / 1024 / 1024;
  }
}

Timers::Timers() : Timers(Verbosity::regular) {}
Timers::Timers(const Verbosity &verbosity)
    : verbosity(verbosity), keep_finished_timers(verbosity >= Verbosity::debug)
{}
Timers::Timers(const Environment &env, const Verbosity &verbosity)
    : Timers(verbosity)
{
  if(env.comm_shared_mem.Rank() != 0)
    {
      // Print info only from the first rank on a node
      this->verbosity = Verbosity::none;
      if(env.num_nodes() != 1)
        node_debug_prefix = El::BuildString("node=", env.node_index(), " ");
    }
}

Timers::~Timers() noexcept
{
  try
    {
      if(verbosity >= Verbosity::debug)
        print_max_mem_used();
    }
  catch(...)
    {
      // destructors should never throw exceptions
    }
}

Timers::Named_Timers::iterator Timers::add_and_start(const std::string &name)
{
  std::string full_name = prefix + name;

  process_meminfo(full_name);

  return named_timers.emplace(named_timers.end(), full_name, Timer());
}
void Timers::write_profile(const std::filesystem::path &path) const
{
  if(path.has_parent_path())
    std::filesystem::create_directories(path.parent_path());
  std::ofstream f(path);

  f << "{" << '\n';
  for(auto it(named_timers.begin()); it != named_timers.end();)
    {
      f << "    {\"" << it->first << "\", " << it->second << "}";
      ++it;
      if(it != named_timers.end())
        {
          f << ",";
        }
      f << '\n';
    }
  f << "}" << '\n';

  ASSERT(f.good(), "Error when writing to: ", path);
}
void Timers::print_max_mem_used() const
{
  if(max_mem_used > 0 && !max_mem_used_name.empty())
    {
      El::Output(node_debug_prefix, "max MemUsed: ", to_GB(max_mem_used),
                 " GB at \"", max_mem_used_name, "\"");
    }
}

// For --verbosity=trace:
// Print memory usage for the current node (from the first rank).
// If we cannot parse /proc/meminfo, then simply print timer name.
//
// In addition, for --verbosity=debug we update max MemUsed (will be printed in the end).
void Timers::process_meminfo(const std::string &name)
{
  // Do not collect memory info for lower verbosity levels, it will not be used
  if(verbosity < Verbosity::debug)
    return;

  const auto prefix = El::BuildString(node_debug_prefix, "start ", name, " ");
  if(!can_read_meminfo)
    {
      // Print "start timer" without MemUsed
      if(verbosity >= Verbosity::trace)
        El::Output(prefix);
      return;
    }

  // Read /proc/meminfo
  constexpr bool print_error_msg = true;
  const auto meminfo
    = Proc_Meminfo::try_read(can_read_meminfo, print_error_msg);

  // Update max MemUsed info, it will be printed for --verbosity=debug
  if(meminfo.mem_used() > max_mem_used)
    {
      max_mem_used = meminfo.mem_used();
      max_mem_used_name = name;
    }

  // Print "start timer", only for --verbosity=trace
  if(verbosity < Verbosity::trace)
    return;

  // MemTotal is constant, thus we print it only once, when adding first timer
  if(named_timers.empty())
    {
      El::Output(prefix, "--- MemTotal: ", to_GB(meminfo.mem_total), " GB");
    }

  //Print MemUsed each time
  El::Output(prefix, "--- MemUsed: ", to_GB(meminfo.mem_used()), " GB");
}
