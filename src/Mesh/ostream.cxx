#include "../Mesh.hxx"
#include "../ostream_array.hxx"

std::ostream & operator<<(std::ostream &os, const Mesh &mesh)
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
