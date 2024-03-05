#include "catch2/catch_amalgamated.hpp"
#include "sdpb_util/Shared_Window_Array.hxx"
#include "test_util/diff.hxx"

#include <El.hpp>
#include <vector>

using Test_Util::REQUIRE_Equal::diff;

TEST_CASE("MPI_Shared_Window")
{
  El::mpi::Comm comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &comm.comm);

  SECTION("Shared_Window_Array")
  {
    size_t size = 10;
    Shared_Window_Array<double> array(comm, size);

    for(size_t i = 0; i < size; ++i)
      {
        if(El::mpi::Rank(comm) == i % El::mpi::Size(comm))
          array[i] = i;
      }

    array.Fence();
    for(size_t i = 0; i < size; ++i)
      {
        DIFF(array[i], (double)i);
      }
  }
}
