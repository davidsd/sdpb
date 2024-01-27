#pragma once

#include "assert.hxx"

#include <El.hpp>
#include <boost/noncopyable.hpp>

template <class T> class Shared_Window_Array : boost::noncopyable
{
public:
  MPI_Win win{};
  El::mpi::Comm comm;
  T *data;
  const size_t size;

public:
  Shared_Window_Array() = delete;
  // shared_memory_comm should be created via
  // MPI_Comm_split_type (MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
  // MPI_INFO_NULL, &shared_memory_comm.comm);
  //
  // It ensures that all ranks in the communicator are on the same node
  // and can share memory.
  Shared_Window_Array(El::mpi::Comm shared_memory_comm, size_t size)
      : comm(shared_memory_comm), size(size)
  {
    MPI_Aint local_window_size; // number of bytes allocated by current rank
    int disp_unit = sizeof(T);

    // Allocate all memory in rank=0
    if(El::mpi::Rank(shared_memory_comm) == 0)
      local_window_size = size * disp_unit;
    else
      local_window_size = 0;

    MPI_Win_allocate_shared(local_window_size, disp_unit, MPI_INFO_NULL,
                            shared_memory_comm.comm, &data, &win);
    // Get local pointer to data allocated in rank=0
    MPI_Win_shared_query(win, 0, &local_window_size, &disp_unit, &data);
    ASSERT(local_window_size == size * sizeof(T));
    ASSERT(disp_unit == sizeof(T));
    Fence();
  }

  ~Shared_Window_Array()
  {
    Fence();
    MPI_Win_free(&win);
  }

  void Fence() const { MPI_Win_fence(0, win); }
  T &operator[](size_t index) { return data[index]; }
  const T &operator[](size_t index) const { return data[index]; }
};
