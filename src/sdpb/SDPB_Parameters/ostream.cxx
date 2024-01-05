#include "../SDPB_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const SDPB_Parameters &p)
{
  os << "SDP directory   : " << p.sdp_path << '\n'
     << "out directory   : " << p.out_directory << '\n'
     << "\nParameters:\n"
     << p.solver
     << "noFinalCheckpoint            = " << p.no_final_checkpoint << '\n'
     << "writeSolution                = " << p.write_solution << '\n'
     << "procGranularity              = " << p.proc_granularity << '\n'
     << "verbosity                    = " << static_cast<int>(p.verbosity)
     << '\n';
  return os;
}
