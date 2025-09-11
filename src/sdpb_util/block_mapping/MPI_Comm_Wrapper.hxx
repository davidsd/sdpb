#pragma once

#include <El.hpp>

struct MPI_Comm_Wrapper
{
  El::mpi::Comm value;
  MPI_Comm_Wrapper() = default;
  // Allow move
  MPI_Comm_Wrapper(MPI_Comm_Wrapper &&other) noexcept = default;
  MPI_Comm_Wrapper &operator=(MPI_Comm_Wrapper &&other) noexcept = default;
  // Prohibit copy
  MPI_Comm_Wrapper(const MPI_Comm_Wrapper &) = delete;
  MPI_Comm_Wrapper &operator=(const MPI_Comm_Wrapper &) = delete;
  ~MPI_Comm_Wrapper()
  {
    if(value != El::mpi::COMM_WORLD)
      {
        El::mpi::Free(value);
      }
  }
};
