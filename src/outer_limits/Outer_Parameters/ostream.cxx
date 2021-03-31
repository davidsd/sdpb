#include "../Outer_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const Outer_Parameters &p)
{
  os << "SDP file   : " << p.sdp_path << '\n'
     << "out directory   : " << p.out_directory << '\n'
     << "\nParameters:\n"
     << p.solver << "noFinalCheckpoint            = " << p.no_final_checkpoint
     << '\n'
     << "writeSolution                = " << p.write_solution << '\n'
     << "verbosity                    = " << static_cast<int>(p.verbosity)
     << '\n';
  return os;
}
