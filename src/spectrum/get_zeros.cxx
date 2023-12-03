#include "sdpb_util/Mesh.hxx"

namespace
{
  void
  maybe_add_point(const El::BigFloat &x_minus, const El::BigFloat &x_bar,
                  const El::BigFloat &x_plus, const El::BigFloat &f_x_minus,
                  const El::BigFloat &f_x_bar, const El::BigFloat &f_x_plus,
                  const El::BigFloat &threshold,
                  std::vector<El::BigFloat> &points)
  {
    const El::BigFloat dx(x_plus - x_minus);
    const El::BigFloat a(f_x_bar), b((f_x_plus - f_x_minus) / dx),
      c((f_x_plus - 2 * f_x_bar + f_x_minus) / (dx * dx / 4));

    // Special case for x=zero, since there the function may go through zero
    if(x_minus == El::BigFloat(0))
      {
      // TODO: This checks if f(0) <= threshold * f'(0).  This has
      // different units from the regular check below.
        
        if(f_x_minus <= threshold * b)
          {
            points.push_back(x_minus);
          }
      }
    else if(c > 0)
      {
        const El::BigFloat x_min(-b / c + x_bar);
        const El::BigFloat f_x_min(a - b * b / (2 * c));
        // This checks if x_min is inside this submesh and
        // f(x_min) < threshold * f''(x_min)
        if(x_min >= x_minus && x_min <= x_plus && f_x_min < c * threshold)
          {
            points.push_back(x_min);
          }
      }
  }
}

std::vector<El::BigFloat>
get_zeros(const Mesh &mesh, const El::BigFloat &threshold)
{
  std::vector<El::BigFloat> result;
  if(mesh.lower)
    {
      for(auto &point : get_zeros(*(mesh.lower), threshold))
        {
          result.push_back(point);
        }
    }
  else
    {
      maybe_add_point(mesh.x[0], mesh.x[1], mesh.x[2], mesh.f[0], mesh.f[1],
                      mesh.f[2], threshold, result);
    }

  if(mesh.upper)
    {
      for(auto &point : get_zeros(*(mesh.upper), threshold))
        {
          result.push_back(point);
        }
    }
  else
    {
      maybe_add_point(mesh.x[2], mesh.x[3], mesh.x[4], mesh.f[2], mesh.f[3],
                      mesh.f[4], threshold, result);
    }
  return result;
}
