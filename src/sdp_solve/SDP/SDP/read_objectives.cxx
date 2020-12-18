#include "../../SDP.hxx"

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void read_objectives(const boost::filesystem::path &sdp_directory,
                     const El::Grid &grid, El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b)
{
  boost::filesystem::ifstream objectives_stream(sdp_directory
                                                / "objectives.json");
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
          dual_objective_b.SetLocal(
            row, 0,
            El::BigFloat(b[global_row].GetString()));
        }
    }
}
