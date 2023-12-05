#include "Mesh.hxx"
#include "sdpb_util/ostream/ostream_array.hxx"

namespace
{
  bool need_refine(const El::BigFloat &f_m, const El::BigFloat &f_x_bar,
                   const El::BigFloat &f_p, const El::BigFloat &mesh_threshold,
                   const El::BigFloat &block_epsilon)
  {
    // TODO: This gets wonky around f()=0
    const El::BigFloat f_bar((f_m + f_p) / 2);
    const El::BigFloat diff(Abs((f_bar - f_x_bar)));
    return diff > mesh_threshold * (Abs(f_bar) + Abs(f_x_bar))
           && diff > block_epsilon;
  }
}

Mesh::Mesh(const El::BigFloat &x_0, const El::BigFloat &x_2,
           const El::BigFloat &x_4, const El::BigFloat &f_0,
           const El::BigFloat &f_2, const El::BigFloat &f_4,
           const std::function<El::BigFloat(const El::BigFloat &x)> &fn,
           const El::BigFloat &mesh_threshold,
           const El::BigFloat &block_epsilon)
    : x({x_0, (x_0 + x_2) / 2, x_2, (x_2 + x_4) / 2, x_4}),
      f({f_0, fn(x[1]), f_2, fn(x[3]), f_4})
{
  // Stop refining if we can not really resolve differences in the coordinates
  // TODO: Scale this by max_delta
  if(Abs(x[0] - x[1]) < Sqrt(El::limits::Epsilon<El::BigFloat>()))
    {
      return;
    }
  if(need_refine(f[0], f[1], f[2], mesh_threshold, block_epsilon))
    {
      lower = std::make_unique<Mesh>(x[0], x[1], x[2], f[0], f[1], f[2], fn,
                                     mesh_threshold, block_epsilon);
    }
  if(need_refine(f[2], f[3], f[4], mesh_threshold, block_epsilon))
    {
      upper = std::make_unique<Mesh>(x[2], x[3], x[4], f[2], f[3], f[4], fn,
                                     mesh_threshold, block_epsilon);
    }
}

std::ostream &operator<<(std::ostream &os, const Mesh &mesh)
{
  os << "{\n  \"x\": " << mesh.x << ",\n  \"f\": " << mesh.f;
  if(mesh.lower)
    {
      os << ",\n    \"lower\": " << *(mesh.lower);
    }
  if(mesh.upper)
    {
      os << ",\n    \"upper\": " << *(mesh.upper);
    }
  os << "}";
  return os;
}
