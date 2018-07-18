#include "read_vector.hxx"

#include <El.hpp>
#include <boost/filesystem.hpp>

void read_objectives(const boost::filesystem::path &sdp_directory,
                     El::BigFloat &objective_const,
                     El::DistMatrix<El::BigFloat> &dual_objective_b)
{
  boost::filesystem::ifstream objectives_stream(sdp_directory / "objectives");
  objectives_stream >> objective_const;

  std::vector<El::BigFloat> temp;
  read_vector(objectives_stream, temp);
  dual_objective_b.Resize(temp.size(), 1);
  if(dual_objective_b.GlobalCol(0) == 0)
    {
      size_t local_height(dual_objective_b.LocalHeight());
      for(size_t row = 0; row < local_height; ++row)
        {
          size_t global_row(dual_objective_b.GlobalRow(row));
          dual_objective_b.SetLocal(row, 0, temp[global_row]);
        }
    }
}
