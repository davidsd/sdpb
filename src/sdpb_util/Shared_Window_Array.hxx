#pragma once

#include "assert.hxx"

#include <El.hpp>
#include <memory>

template <class T> class Abstract_Shared_Window_Array
{
public:
  virtual ~Abstract_Shared_Window_Array() = default;

  [[nodiscard]] virtual T *data() = 0;
  [[nodiscard]] virtual size_t size() const = 0;

  virtual void Fence() const = 0;
  T &operator[](size_t index) { return data()[index]; }
  const T &operator[](size_t index) const { return data()[index]; }
};

template <class T>
class Shared_Window_Array final : public Abstract_Shared_Window_Array<T>
{
  MPI_Win _win{};
  El::mpi::Comm _comm;
  T *_data;
  size_t _size = 0;

public:
  Shared_Window_Array(const Shared_Window_Array &other) = delete;
  Shared_Window_Array(Shared_Window_Array &&other) noexcept = default;
  Shared_Window_Array &operator=(const Shared_Window_Array &other) = delete;
  Shared_Window_Array &operator=(Shared_Window_Array &&other) noexcept
    = default;

  Shared_Window_Array() = default;
  // shared_memory_comm should be created via
  // MPI_Comm_split_type (MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
  // MPI_INFO_NULL, &shared_memory_comm.comm);
  //
  // It ensures that all ranks in the communicator are on the same node
  // and can share memory.
  Shared_Window_Array(El::mpi::Comm shared_memory_comm, size_t size)
      : _comm(shared_memory_comm), _size(size)
  {
    MPI_Aint local_window_size; // number of bytes allocated by current rank
    int disp_unit = sizeof(T);

    // Allocate all memory in rank=0
    if(El::mpi::Rank(shared_memory_comm) == 0)
      local_window_size = size * disp_unit;
    else
      local_window_size = 0;

    MPI_Win_allocate_shared(local_window_size, disp_unit, MPI_INFO_NULL,
                            shared_memory_comm.comm, &_data, &_win);
    // Get local pointer to data allocated in rank=0
    MPI_Win_shared_query(_win, 0, &local_window_size, &disp_unit, &_data);
    ASSERT_EQUAL(local_window_size, size * sizeof(T));
    ASSERT_EQUAL(disp_unit, sizeof(T));
    Shared_Window_Array::Fence();
  }

  ~Shared_Window_Array() override
  {
    // If one rank throws an exception and another doesn't,
    // The first rank will wait (potentially forever)
    // at the fence instead of exiting.
    // To prevent it, we disable synchronization if an exception has been thrown.
    // NB: if exception is caught after that and program continues working,
    // it will probably hang on the next synchronization point!
    if(std::uncaught_exceptions() == 0 || _comm.Size() == 1)
      {
        Shared_Window_Array::Fence();
        MPI_Win_free(&_win);
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

  T *data() override { return _data; }
  size_t size() const override { return _size; }
  void Fence() const override { MPI_Win_fence(0, _win); }
  [[nodiscard]] El::mpi::Comm Comm() const { return _comm; }
};

template <class T>
class Shared_Window_Array_View final : public Abstract_Shared_Window_Array<T>
{
  std::shared_ptr<Shared_Window_Array<T>> _array;
  const size_t _offset;
  const size_t _size;

public:
  // View over existing shared memory window
  Shared_Window_Array_View(std::shared_ptr<Shared_Window_Array<T>> array,
                           const size_t offset, const size_t size)
      : _array(array), _offset(offset), _size(size)
  {
    ASSERT(
      array->size() >= offset + size,
      "Shared_Window_Array_View out of bounds:", DEBUG_STRING(array->size()),
      DEBUG_STRING(offset), DEBUG_STRING(size));
  }
  // View through another view
  Shared_Window_Array_View(Shared_Window_Array_View &array_view,
                           const size_t offset, const size_t size)
      : Shared_Window_Array_View(array_view._array,
                                 array_view._offset + offset, size)
  {
    ASSERT(size <= array_view.size(),
           "Requested Shared_Window_Array_View size=", size,
           " is larger than the parent view size=", array_view.size());
  }
  // Create new shared memory window and return a view over it
  Shared_Window_Array_View(El::mpi::Comm shared_memory_comm, size_t size)
      : Shared_Window_Array_View(
          std::make_shared<Shared_Window_Array<T>>(shared_memory_comm, size),
          0, size)
  {}
  T *data() override { return _array->data() + _offset; }
  size_t size() const override { return _size; }
  [[nodiscard]] El::mpi::Comm Comm() const { return _array->Comm(); }
  void Fence() const override { _array->Fence(); }
};
