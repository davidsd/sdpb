#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

namespace fs = std::filesystem;

void write_profiling(const fs::path &spectrum_output_path, Timers &timers)
{
  Scoped_Timer write_profiling_timer(timers, "write_profiling");

  fs::path profiling_dir = spectrum_output_path;
  profiling_dir += ".profiling";
  if(El::mpi::Rank() == 0)
    {
      std::error_code ec;
      fs::create_directories(profiling_dir, ec);
      if(ec)
        PRINT_WARNING("Failed to create directory ", profiling_dir, ": ",
                      ec.message());
    }
  try
    {
      // Wait until directory is created
      El::mpi::Barrier();
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
