#include "../Outer_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const Outer_Parameters &p)
{
  os << "functions file  : " << p.functions_path << '\n'
     << "out directory   : " << p.output_path << '\n'
     << "\nParameters:\n"
     << "dualityGapReduction          = " << p.duality_gap_reduction << '\n'
     << "meshThreshold                = " << p.mesh_threshold << '\n'
     << p.solver
     << "writeSolution                = " << p.write_solution << '\n'
     << "verbosity                    = " << static_cast<int>(p.verbosity)
     << '\n';
  return os;
}
