#include "Proc_Meminfo.hxx"

#include <El.hpp>

#include <fstream>
#include <sstream>

Proc_Meminfo::Proc_Meminfo(size_t mem_total, size_t mem_available)
    : mem_total(mem_total), mem_available(mem_available)
{}

size_t Proc_Meminfo::mem_used() const
{
  return mem_total - mem_available;
}

Proc_Meminfo
Proc_Meminfo::try_read(bool &result, bool print_error_msg) noexcept
{
  try
    {
      const auto meminfo = read();
      result = true;
      return meminfo;
    }
  catch(std::exception &e)
    {
      result = false;
      if(print_error_msg)
        El::Output("Failed to parse /proc/meminfo: ", e.what());
      return {0, 0};
    }
  catch(...)
    {
      result = false;
      if(print_error_msg)
        El::Output("Failed to parse /proc/meminfo");
      return {0, 0};
    }
}

Proc_Meminfo Proc_Meminfo::read() noexcept(false)
{
  auto proc_meminfo_path = "/proc/meminfo";
  std::ifstream meminfo_file(proc_meminfo_path);

  if(!meminfo_file.good())
    El::RuntimeError("Cannot open ", proc_meminfo_path);

  const char *mem_total_prefix = "MemTotal:";
  const char *mem_available_prefix = "MemAvailable:";
  size_t memTotalKB = 0;
  size_t memAvailableKB = 0;
  std::string line;
  while(std::getline(meminfo_file, line))
    {
      std::istringstream iss(line);
      std::string name;
      if(iss >> name)
        {
          if(name != mem_total_prefix && name != mem_available_prefix)
            continue;
          size_t size;
          std::string kB;
          if(iss >> size >> kB)
            {
              if(kB != "kB" && kB != "KB")
                {
                  El::RuntimeError(
                    proc_meminfo_path,
                    ": expected \"kB\" at the end of line: ", line);
                }
              if(name == mem_total_prefix)
                memTotalKB = size;
              else if(name == mem_available_prefix)
                memAvailableKB = size;
              if(memTotalKB > 0 && memAvailableKB > 0)
                break;
              continue;
            }
        }
      El::RuntimeError(proc_meminfo_path, ": cannot parse line: ", line);
    }

  if(memTotalKB == 0)
    {
      El::RuntimeError(proc_meminfo_path, ": ", mem_total_prefix,
                       " not found");
    }
  if(memAvailableKB == 0)
    {
      El::RuntimeError(proc_meminfo_path, ": ", mem_available_prefix,
                       " not found");
    }

  return {memTotalKB * 1024, memAvailableKB * 1024};
}
