#include "../Dynamical_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const Dynamical_Parameters &p)
{
  os << "centering SDP directory   : " << p.sdp_path
     << '\n'
     //<< "shifted SDP directory  : "    << p.new_sdp_path << '\n'
     << "out directory   : " << p.out_directory << '\n'
     << "\nParameters:\n"
     << p.solver << "noFinalCheckpoint            = " << p.no_final_checkpoint
     << '\n'
     << "writeSolution                = " << p.write_solution << '\n'
     << "procGranularity              = " << p.proc_granularity << '\n'
     << "verbosity                    = " << static_cast<int>(p.verbosity)
     << '\n';

  std::cout << "search_direction             = ";
  for(auto i : p.solver.search_direction)
    std::cout << i << ' ';
  std::cout << "\n\n";

  return os;
}
