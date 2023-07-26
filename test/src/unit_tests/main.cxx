#include <catch2/catch_amalgamated.hpp>
#include <El.hpp>
#include <iostream>
#include <sstream>

#ifndef CATCH_AMALGAMATED_CUSTOM_MAIN
#error "To override main, pass '-D CATCH_AMALGAMATED_CUSTOM_MAIN' to compiler"
#endif

int main(int argc, char *argv[])
{
  El::Environment env(argc, argv);
  El::gmp::SetPrecision(128);

  int rank = El::mpi::Rank();

  // redirect cout to string stream
  std::stringstream ss;
  auto cout_buf = std::cout.rdbuf(ss.rdbuf());

  int result = Catch::Session().run(argc, argv);

  // restore cout buffer
  std::cout.rdbuf(cout_buf);

  // Print output.
  // Always print from master rank,
  // print from other ranks only if they failed.
  for(int i = 0; i < El::mpi::Size(); ++i)
    {
      // Ensure rank printing order with barrier.
      // TODO sometimes output is still shuffled.
      // Better write to file and then print all from rank 0?
      El::mpi::Barrier();
      if(i == rank)
        {
          if(rank == 0 || result != 0)
            {
              El::Output("MPI Rank ", rank);
              El::Output(ss.str());
            }
        }
    }
  return result;
}
