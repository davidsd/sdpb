#pragma once

#include <El.hpp>

struct MPI_Comm_Wrapper
{
  El::mpi::Comm value;
  MPI_Comm_Wrapper() = default;
  MPI_Comm_Wrapper(const MPI_Comm_Wrapper &) = delete;
  void operator=(const MPI_Comm_Wrapper &) = delete;
  ~MPI_Comm_Wrapper()
  {
    if(value != El::mpi::COMM_WORLD)
      {
        El::mpi::Free(value);
      }
  }
};

namespace std
{
  inline void swap(MPI_Comm_Wrapper &a, MPI_Comm_Wrapper &b) noexcept
  {
    swap(a.value, b.value);
  }
}
