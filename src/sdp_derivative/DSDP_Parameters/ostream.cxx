#include "../DSDP_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const DSDP_Parameters &p)
{
  os << std::boolalpha
     << "sdp                          = " << p.sdp_path << '\n'
     << "newSdp                       = " << p.new_sdp_path << '\n'
     << "procsPerNode                 = " << p.procs_per_node << '\n'
     << "procsGranularity             = " << p.proc_granularity << '\n'
     << "precision(actual)            = " << p.precision << "("
     << mpf_get_default_prec() << ")" << '\n'

     << "outDir                       = " << p.out_directory << '\n'
     << "solutionDir                  = " << p.solution_dir << '\n'
     << "verbosity                    = " << static_cast<int>(p.verbosity)
     << '\n';
  return os;
}
