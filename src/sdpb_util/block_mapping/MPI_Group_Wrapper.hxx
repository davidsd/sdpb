#pragma once

#include <El.hpp>

struct MPI_Group_Wrapper
{
  El::mpi::Group value;

  MPI_Group_Wrapper() = default;
  // Allow move
  MPI_Group_Wrapper(MPI_Group_Wrapper &&other) noexcept = default;
  MPI_Group_Wrapper &operator=(MPI_Group_Wrapper &&other) noexcept = default;
  // Prohibit copy
  MPI_Group_Wrapper(const MPI_Group_Wrapper &) = delete;
  MPI_Group_Wrapper &operator=(const MPI_Group_Wrapper &) = delete;

  ~MPI_Group_Wrapper()
  {
    if(value != El::mpi::GROUP_NULL)
      {
        El::mpi::Free(value);
      }
  }
};
