#include "../Approx_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const Approx_Parameters &p)
{
  os << std::boolalpha
     << "sdp                          = " << p.sdp_path << '\n'
     << "newSdp                       = " << p.new_sdp_path << '\n'
     << "procsGranularity             = " << p.proc_granularity << '\n'
     << "precision(actual)            = " << p.precision << "("
     << mpf_get_default_prec() << ")" << '\n'
     << "solutionDir                  = " << p.solution_dir << '\n'
     << "writeSolverState             = " << p.write_solver_state << '\n'
     << '\n';
  return os;
}
