#include "../sdp_solve.hxx"

El::Matrix<El::BigFloat>
read_y(const boost::filesystem::path &solution_dir, const size_t &y_height)
{
  El::DistMatrix<El::BigFloat> y_dist(y_height, 1);
  read_text_block(y_dist, solution_dir / "y.txt");
  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> y_star(y_dist);

  El::Matrix<El::BigFloat> result(y_star.LockedMatrix());
  return result;
}
