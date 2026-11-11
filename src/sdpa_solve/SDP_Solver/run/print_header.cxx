#include "sdpa_solve/SDP_Solver.hxx"

#include <iostream>

namespace Sdpb::Sdpa
{
  void print_header(const Verbosity &verbosity)
  {
    if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
      {
        std::cout
          << "\n"
          << "          time    mu     P-obj       D-obj      gap     "
             "    P-err       D-err      P-step   D-step   beta\n"
          << "--------------------------------------------------------"
             "-------------------------------------------------------------\n";
      }
  }
}
