#include "../../Dynamical_Solver.hxx"

#include <iostream>

void print_header_dynamical(const Verbosity &verbosity)
{
  if(verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      std::cout << "\n"
                << "          time    mu     P-obj       D-obj      gap     "
                   "    P-err       p-err       D-err      P-step   D-step   beta   ext-step-size\n"
                << "--------------------------------------------------------"
                   "-----------------------------------------------------------------------------\n";
    }
}

void print_header_dynamical_new(const Verbosity &verbosity)
{
	if (verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
	{
		std::cout << "\n"
			<< "          time    mu     P-obj       D-obj      gap     "
			"    P-err       p-err        D-err        R-err      P-step   D-step   beta   ext-step-size\n"
			<< "--------------------------------------------------------"
			"----------------------------------------------------------------------------------\n";
	}
}

