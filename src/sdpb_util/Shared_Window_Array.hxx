#pragma once

#include "assert.hxx"

#include <El.hpp>

template <class T> class Shared_Window_Array
{
public:
  Shared_Window_Array(const Shared_Window_Array &other) = delete;
  Shared_Window_Array(Shared_Window_Array &&other) noexcept = default;
  Shared_Window_Array &operator=(const Shared_Window_Array &other) = delete;
  Shared_Window_Array &operator=(Shared_Window_Array &&other) noexcept
    = default;

  MPI_Win win{};
  El::mpi::Comm comm;
  T *data;
  size_t size = 0;

public:
  Shared_Window_Array() = default;
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
    ASSERT_EQUAL(local_window_size, size * sizeof(T));
    ASSERT_EQUAL(disp_unit, sizeof(T));
    Fence();
  }

  ~Shared_Window_Array()
  {
    // If one rank throws an exception and another doesn't,
    // The first rank will wait (potentially forever)
    // at the fence instead of exiting.
    // To prevent it, we disable synchronization if an exception has been thrown.
    // NB: if exception is caught after that and program continues working,
    // it will probably hang on the next synchronization point!
    if(std::uncaught_exceptions() == 0 || comm.Size() == 1)
      {
        Fence();
        MPI_Win_free(&win);
      }
    else
      {
        try
          {
            PRINT_WARNING(
              "~Shared_Window_Array() called during stack unwinding on rank=",
              El::mpi::Rank(),
              ". The program should exit after that, otherwise MPI will not "
              "work correctly.");
          }
        catch(...)
          {}
      }
  }

  void Fence() const { MPI_Win_fence(0, win); }
  T &operator[](size_t index) { return data[index]; }
  const T &operator[](size_t index) const { return data[index]; }
};
