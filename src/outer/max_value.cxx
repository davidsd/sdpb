#include "Mesh.hxx"

El::BigFloat max_value(const Mesh &mesh)
{
  El::BigFloat result(0);
  for(auto &f: mesh.f)
    {
      result=Max(result,Abs(f));
    }
  if(mesh.lower)
    {
      result=Max(result,max_value(*(mesh.lower)));
    }
  if(mesh.upper)
    {
      result=Max(result,max_value(*(mesh.upper)));
    }
  return result;
}

