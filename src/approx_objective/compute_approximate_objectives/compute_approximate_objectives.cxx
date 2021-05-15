#include "../../sdp_solve.hxx"
#include "../../sdp_read.hxx"

#include <boost/filesystem.hpp>

void Axpy(const El::BigFloat &alpha, const SDP &new_sdp, SDP &delta_sdp);

El::BigFloat compute_approximate_objective(
  const Block_Info &block_info, const SDP &sdp, const SDP &d_sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q);

std::vector<El::BigFloat> compute_approximate_objectives(
  const Block_Info &block_info, const El::Grid &grid, const SDP &sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q,
  const boost::filesystem::path &input_path)
{
  std::vector<El::BigFloat> result;
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

      result.push_back(compute_approximate_objective(
        block_info, sdp, d_sdp, x, y, schur_complement_cholesky,
        schur_off_diagonal, Q));
    }
  return result;
}
