#pragma once

#include <El.hpp>

struct MPI_Group_Wrapper
{
  El::mpi::Group value;
  MPI_Group_Wrapper() = default;
  MPI_Group_Wrapper(const MPI_Group_Wrapper &) = delete;
  void operator=(const MPI_Group_Wrapper &) = delete;
  ~MPI_Group_Wrapper()
  {
    if(value != El::mpi::GROUP_NULL)
      {
        El::mpi::Free(value);
      }
  }
};

namespace std
{
  inline void swap(MPI_Group_Wrapper &a, MPI_Group_Wrapper &b) noexcept
  {
    swap(a.value, b.value);
  }
}
