#include "../Dynamical_Solver_Parameters.hxx"

std::ostream &operator<<(std::ostream &os, const Dynamical_Solver_Parameters &p)
{
  os << "stepSizeAlpha                = " << p.alpha << '\n'
     << "updateSdpThresholdMax        = " << p.update_sdp_threshold_max << '\n'
     << "updateSdpThresholdMin        = " << p.update_sdp_threshold_min << '\n'
     << "numExternalParams            = " << p.n_external_parameters << '\n' 
     << "shiftedSDPsDir               = " << p.new_sdp_path << '\n'
     << "findBoundary                 = " << p.find_boundary << '\n'
     << "findBoundaryObjThreshold     = " << p.find_boundary_obj_threshold << '\n'
     << "muDirectionMode              = " << p.mu_last_direction << '\n'
     << p.solver_parameters 
     << '\n';
 
  return os;
}
