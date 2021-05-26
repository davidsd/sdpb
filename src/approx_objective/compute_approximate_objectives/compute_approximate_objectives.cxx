#include "../Approx_Objective.hxx"
#include "../../sdp_solve.hxx"
#include "../../sdp_read.hxx"

#include <boost/filesystem.hpp>

void Axpy(const El::BigFloat &alpha, const SDP &new_sdp, SDP &delta_sdp);

Approx_Objective compute_approximate_objective(
  const Block_Info &block_info, const SDP &sdp, const SDP &d_sdp,
  const El::BigFloat &new_objective_const, const Block_Vector &x,
  const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q);

std::vector<std::pair<std::string, Approx_Objective>>
compute_approximate_objectives(
  const Block_Info &block_info, const El::Grid &grid, const SDP &sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q,
  const boost::filesystem::path &input_path)
{
  std::vector<std::pair<std::string, Approx_Objective>> result;
  if(input_path.extension() == ".nsv")
    {
      for(auto &filename : read_file_list(input_path))
        {
          for(auto &objective : compute_approximate_objectives(
                block_info, grid, sdp, x, y, schur_complement_cholesky,
                schur_off_diagonal, Q, filename))
            {
              result.push_back(objective);
            }
        }
    }
  else
    {
      SDP new_sdp(input_path, block_info, grid), d_sdp(new_sdp);
      Axpy(El::BigFloat(-1), sdp, d_sdp);

      result.emplace_back(input_path.string(),
                          compute_approximate_objective(
                            block_info, sdp, d_sdp, new_sdp.objective_const, x,
                            y, schur_complement_cholesky, schur_off_diagonal,
                            Q));
    }
  return result;
}
