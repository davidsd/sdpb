#include "../../../read_vector.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void read_objectives(const boost::filesystem::path &sdp_directory,
                     const El::Grid &grid, El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objectives_b)
{
  const boost::filesystem::path objectives_path(sdp_directory / "objectives");
  boost::filesystem::ifstream objectives_stream(objectives_path);
  if(!objectives_stream.good())
    {
      throw std::runtime_error("Could not open '" + objectives_path.string()
                               + "'");
    }
  objectives_stream >> objective_const;
  if(!objectives_stream.good())
    {
      throw std::runtime_error("Corrupted file: " + objectives_path.string());
    }

  std::vector<El::BigFloat> temp;
  read_vector(objectives_stream, temp);
  dual_objectives_b.SetGrid(grid);
  dual_objectives_b.Resize(temp.size(), 1);
  if(dual_objectives_b.GlobalCol(0) == 0)
    {
      size_t local_height(dual_objectives_b.LocalHeight());
      for(size_t row = 0; row < local_height; ++row)
        {
          size_t global_row(dual_objectives_b.GlobalRow(row));
          dual_objectives_b.SetLocal(row, 0, temp[global_row]);
        }
    }
}
