#include "../../SDP.hxx"

#include <archive_reader.hpp>
#include <archive_exception.hpp>

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

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

void read_objectives(const boost::filesystem::path &sdp_path,
                     const El::Grid &grid, El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b)
{
  const std::string objectives_name("objectives.json");
  if(boost::filesystem::is_regular_file(sdp_path))
    {
      // TODO: This is going to reopen the zip file many, many
      // times.
      boost::filesystem::ifstream fs(sdp_path);
      ns_archive::reader reader
        (ns_archive::reader::make_reader<ns_archive::ns_reader::format::_ALL,
                                          ns_archive::ns_reader::filter::_ALL>(
                                                                               fs, 10240));

      for(auto entry : reader)
        {
          if(entry->get_header_value_pathname() == objectives_name)
            {
              read_objectives_stream(grid, entry->get_stream(),
                                     objective_const, dual_objective_b);
            }
        }
    }
  else
    {
      boost::filesystem::ifstream objectives_stream(sdp_path
                                                    / objectives_name);
      read_objectives_stream(grid, objectives_stream, objective_const,
                             dual_objective_b);
    }
}
