#include "../Dynamic_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const Dynamic_Parameters &p)
{
  os << "centering SDP directory   : " << p.sdp_path << '\n'
     << "surronding SDP directory  : " << p.new_sdp_path << '\n'
     << "out directory   : " << p.out_directory << '\n'
     << "\nParameters:\n"
     << p.solver
     << "noFinalCheckpoint            = " << p.no_final_checkpoint << '\n'
     << "writeSolution                = " << p.write_solution << '\n'
     << "procsPerNode                 = " << p.procs_per_node << '\n'
     << "procGranularity              = " << p.proc_granularity << '\n'
     << "verbosity                    = " << static_cast<int>(p.verbosity) << '\n'
     << "StepSizeAlpha                = " << p.alpha << '\n'
     << "NumberExternalParameters     = " << p.N_external_parameters
     << '\n';
  return os;
}
