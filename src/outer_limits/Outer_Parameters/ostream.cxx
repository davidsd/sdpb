#include "../Outer_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const Outer_Parameters &p)
{
  os << "SDP file   : " << p.sdp_path << '\n'
     << "out directory   : " << p.output_path << '\n'
     << "\nParameters:\n"
     << "dualityGapReduction          = " << p.duality_gap_reduction << '\n'
     << p.solver
     << "writeSolution                = " << p.write_solution << '\n'
     << "verbosity                    = " << static_cast<int>(p.verbosity)
     << '\n';
  return os;
}