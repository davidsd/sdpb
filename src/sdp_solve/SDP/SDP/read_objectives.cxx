#include "sdp_solve/SDP.hxx"
#include "sdp_solve/Archive_Reader.hxx"

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <filesystem>

namespace fs = std::filesystem;

namespace
{
  void
  read_objectives_stream(const El::Grid &grid, std::istream &objectives_stream,
                         El::BigFloat &objective_const,
                         El::DistMatrix<El::BigFloat> &dual_objective_b)
  {
    rapidjson::IStreamWrapper wrapper(objectives_stream);
    rapidjson::Document d;
    d.ParseStream(wrapper);
    objective_const = El::BigFloat(d["constant"].GetString());
    auto b(d["b"].GetArray());

    dual_objective_b.SetGrid(grid);
    dual_objective_b.Resize(b.Size(), 1);
    if(dual_objective_b.GlobalCol(0) == 0)
      {
        size_t local_height(dual_objective_b.LocalHeight());
        for(size_t row = 0; row < local_height; ++row)
          {
            size_t global_row(dual_objective_b.GlobalRow(row));
            dual_objective_b.SetLocal(row, 0,
                                      El::BigFloat(b[global_row].GetString()));
          }
      }
  }
}

void read_objectives(const fs::path &sdp_path, const El::Grid &grid,
                     El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b,
                     Timers &timers)
{
  Scoped_Timer timer(timers, "read_objectives");
  if(!fs::exists(sdp_path))
    {
      El::RuntimeError("SDP path '" + sdp_path.string() + "' does not exist");
    }

  const std::string objectives_name("objectives.json");
  if(fs::is_regular_file(sdp_path))
    {
      // TODO: This is going to reopen the zip file many, many
      // times.
      Archive_Reader reader(sdp_path);
      while(reader.next_entry())
        {
          if(objectives_name == archive_entry_pathname(reader.entry_ptr))
            {
              std::istream stream(&reader);
              read_objectives_stream(grid, stream,
                                     objective_const, dual_objective_b);
            }
        }
      if(dual_objective_b.Height() == 0)
        {
          throw std::runtime_error("Unable to read objectives from input");
        }
    }
  else
    {
      std::ifstream objectives_stream(sdp_path / objectives_name);
      read_objectives_stream(grid, objectives_stream, objective_const,
                             dual_objective_b);
    }
}
