#include "../Dynamical_Solver_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const Dynamical_Solver_Parameters &p)
{
  os << "stepSizeAlpha                = " << p.alpha << '\n'
     << "updateSdpThresholdMax        = " << p.update_sdp_threshold_max << '\n'
     << "numExternalParams            = " << p.n_external_parameters << '\n' 
     << "shiftedSDPsDir               = " << p.new_sdp_path << '\n'
     << p.solver_parameters 
     << '\n';
 
  return os;
}
