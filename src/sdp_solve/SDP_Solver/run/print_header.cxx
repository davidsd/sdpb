#include "../../SDP_Solver.hxx"

#include <iostream>

void print_header(const Verbosity &verbosity)
{
  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      std::cout << "\n"
                << "          time    mu     P-obj       D-obj      gap     "
                   "    P-err       p-err       D-err      P-step   D-step   beta\n"
                << "--------------------------------------------------------"
                   "-------------------------------------------------------------\n";
    }
}
