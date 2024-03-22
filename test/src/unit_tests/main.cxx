#include "sdpb_util/Environment.hxx"

#include <catch2/catch_amalgamated.hpp>
#include <El.hpp>
#include <iostream>
#include <sstream>
#include "test_util/diff.hxx"

#ifndef CATCH_AMALGAMATED_CUSTOM_MAIN
#error "To override main, pass '-D CATCH_AMALGAMATED_CUSTOM_MAIN' to compiler"
#endif

int main(int argc, char *argv[])
{
  Environment env(argc, argv);
  Environment::set_precision(768);
  Test_Util::REQUIRE_Equal::diff_precision = El::gmp::Precision();

  int rank = El::mpi::Rank();

  int result;
  std::stringstream ss;

  if(rank == 0)
    {
      El::Output("MPI Rank ", rank);
      result = Catch::Session().run(argc, argv);
    }
  else
    {
      // redirect cout to string stream
      auto cout_buf = std::cout.rdbuf(ss.rdbuf());
      // run
      result = Catch::Session().run(argc, argv);
      // restore cout buffer
      std::cout.rdbuf(cout_buf);
    }

  // We've printed output from master rank,
  // Now let's print from all other ranks that failed.
  for(int i = 1; i < El::mpi::Size(); ++i)
    {
      // Ensure rank printing order with barrier.
      // TODO sometimes output is still shuffled.
      // Better write to file and then print all from rank 0?
      El::mpi::Barrier();
      if(i == rank)
        {
          if(result != 0)
            {
              El::Output("MPI Rank ", rank);
              El::Output(ss.str());
            }
        }
    }
  return result;
}
