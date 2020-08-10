#include "set_dual_objective_b.hxx"
#include "../../read_vector.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

void read_objectives(const boost::filesystem::path &sdp_directory,
                     const El::Grid &grid, El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b)
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
  set_dual_objective_b(temp, grid, dual_objective_b);
}
