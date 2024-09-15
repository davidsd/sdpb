#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

namespace fs = std::filesystem;

namespace
{
  fs::path get_profiling_dir(const fs::path &spectrum_output_path)
  {
    fs::path profiling_dir = spectrum_output_path;
    profiling_dir += ".profiling";
    return profiling_dir;
  }
}

void create_profiling_dir(const fs::path &spectrum_output_path)
{
  const auto profiling_dir = get_profiling_dir(spectrum_output_path);
  if(El::mpi::Rank() == 0)
    {
      std::error_code ec;
      fs::create_directories(profiling_dir, ec);
      if(ec)
        PRINT_WARNING("Failed to create directory ", profiling_dir, ": ",
                      ec.message());
    }
  // Wait until directory is created
  El::mpi::Barrier();
}

void write_profiling(const fs::path &spectrum_output_path, Timers &timers)
{
  Scoped_Timer write_profiling_timer(timers, "write_profiling");

  const auto profiling_dir = get_profiling_dir(spectrum_output_path);
  try
    {
      ASSERT(fs::exists(profiling_dir) && fs::is_directory(profiling_dir),
             DEBUG_STRING(profiling_dir));
      timers.write_profile(profiling_dir
                           / ("profiling." + std::to_string(El::mpi::Rank())));
    }
  catch(const std::exception &e)
    {
      PRINT_WARNING("Failed to write profiling data to ", profiling_dir, ": ",
                    e.what());
    }
}
