#include "../Mesh.hxx"

bool is_origin_zero(const Mesh &mesh, const El::BigFloat &threshold)
{
  if(mesh.lower)
    {
      return is_origin_zero(*(mesh.lower), threshold);
    }
  else
    {
      // TODO: This checks if f(0) < threshold * f'(0).  This has
      // different units from the regular check, which is f(x) <
      // threshold * f''(x).
      return (mesh.f[0]
              < threshold * (mesh.f[1] - mesh.f[0]) / (mesh.x[1] - mesh.x[0]));
    }
}
