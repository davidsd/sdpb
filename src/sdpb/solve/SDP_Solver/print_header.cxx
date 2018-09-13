#include "../SDP_Solver.hxx"

#include <iostream>

void SDP_Solver::print_header()
{
  if(El::mpi::Rank() == 0)
    {
      std::cout << "\n          time    mu     P-obj       D-obj      gap     "
                   "    P-err       D-err      P-step   D-step   beta  dim\n";
      std::cout
        << "------------------------------------------------------------------"
           "---------------------------------------------\n";
    }
}
