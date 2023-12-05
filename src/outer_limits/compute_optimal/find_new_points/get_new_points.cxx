#include "sdpb_util/Mesh.hxx"

void maybe_add_point(const El::BigFloat &x_minus, const El::BigFloat &x_bar,
                     const El::BigFloat &x_plus, const El::BigFloat &f_x_minus,
                     const El::BigFloat &f_x_bar, const El::BigFloat &f_x_plus,
                     const El::BigFloat &block_epsilon,
                     std::vector<El::BigFloat> &points)
{
  const El::BigFloat dx(x_plus - x_minus);
  const El::BigFloat a(f_x_bar), b((f_x_plus - f_x_minus) / dx),
    c((f_x_plus - 2 * f_x_bar + f_x_minus) / (dx * dx / 4));

  if(c > 0)
    {
      const El::BigFloat x_min(-b / c + x_bar);
      const El::BigFloat f_x_min(a - b * b / (2 * c)),
        f_bar((f_x_plus + f_x_minus) / 2);
      if(x_min >= x_minus && x_min <= x_plus && f_x_min < Abs(f_x_bar - f_bar)
         && Abs(f_x_min) > block_epsilon)
        {
          points.push_back(x_min);
        }
    }
}

std::vector<El::BigFloat>
get_new_points(const Mesh &mesh, const El::BigFloat &block_epsilon)
{
  std::vector<El::BigFloat> result;
  if(mesh.lower)
    {
      for(auto &point : get_new_points(*(mesh.lower), block_epsilon))
        {
          result.push_back(point);
        }
    }
  else
    {
      maybe_add_point(mesh.x[0], mesh.x[1], mesh.x[2], mesh.f[0], mesh.f[1],
                      mesh.f[2], block_epsilon, result);
    }

  if(mesh.upper)
    {
      for(auto &point : get_new_points(*(mesh.upper), block_epsilon))
        {
          result.push_back(point);
        }
    }
  else
    {
      maybe_add_point(mesh.x[2], mesh.x[3], mesh.x[4], mesh.f[2], mesh.f[3],
                      mesh.f[4], block_epsilon, result);
    }
  return result;
}
