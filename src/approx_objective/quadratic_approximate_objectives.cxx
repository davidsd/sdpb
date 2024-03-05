#include "Approx_Objective.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "pmp_read/pmp_read.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void Axpy(const El::BigFloat &alpha, const SDP &new_sdp, SDP &delta_sdp);

std::vector<std::pair<std::string, Approx_Objective>>
quadratic_approximate_objectives(
  const Block_Info &block_info, const El::Grid &grid, const SDP &sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q, const fs::path &input_path)
{
  std::vector<std::pair<std::string, Approx_Objective>> result;
  if(input_path.extension() == ".nsv")
    {
      for(auto &filename : read_nsv_file_list(input_path))
        {
          for(auto &objective : quadratic_approximate_objectives(
                block_info, grid, sdp, x, y, schur_complement_cholesky,
                schur_off_diagonal, Q, filename))
            {
              result.push_back(objective);
            }
        }
    }
  else
    {
      Timers timers;
      SDP new_sdp(input_path, block_info, grid, timers), d_sdp(new_sdp);
      Axpy(El::BigFloat(-1), sdp, d_sdp);

      result.emplace_back(input_path.string(),
                          Approx_Objective(block_info, sdp, d_sdp, x, y,
                                           schur_complement_cholesky,
                                           schur_off_diagonal, Q));
    }
  return result;
}
