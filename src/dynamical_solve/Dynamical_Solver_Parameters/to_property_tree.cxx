#include "../Dynamical_Solver_Parameters.hxx"

boost::property_tree::ptree to_property_tree(const Dynamical_Solver_Parameters &p)
{
  boost::property_tree::ptree result(to_property_tree(p.solver_parameters));
  result.put("newSdpDirs", p.new_sdp_path.string()); 
  result.put("stepSizeAlpha", p.alpha);
  result.put("numExternalParams", p.n_external_parameters);
  result.put("updateSdpThreshold", p.update_sdp_threshold_max);
 

  return result;
}
