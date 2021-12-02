#include "../Dynamical_Solver_Parameters.hxx"

#include <boost/program_options.hpp>

boost::program_options::options_description Dynamical_Solver_Parameters::options()
{
  boost::program_options::options_description result("Dynamic Solver parameters");
  // We set default parameters using El::BigFloat("1e-10",10)
  // rather than a straight double precision 1e-10 so that results
  // are more reproducible at high precision.  Using double
  // precision defaults results in differences of about 1e-15 in
  // primalObjective after one step.

  result.add_options()("newSdpDirs", 
                              boost::program_options::value<boost::filesystem::path>(&new_sdp_path)->required(),
                              "Directory containing the preprocessed SDP data files around the center SDP in external parameter space.");  
 
  result.add_options()("stepSizeAlpha", 
                              boost::program_options::value<El::BigFloat>(&alpha)->default_value(1),
                              "Step size in the external-parameter space to generate the new SDP data files. "
                              "The default value is set to 1.");
  result.add_options()("numExternalParams",
                              boost::program_options::value<int>(&n_external_parameters)->default_value(0),
                              "The number of external parameters to be varied in each iteration of the dynamical SDP. "
                              "The default value is set to 0.");
  result.add_options()("updateSdpThresholdMax",
                              boost::program_options::value<El::BigFloat>(&update_sdp_threshold_max)->default_value(1),
                              "Take a step in the external parameters, "
                              "that is to regenerator the sdp files if the step size is smaller than the threshold. "
                              "The default value is set to 1.");  
  result.add(solver_parameters.options()); 

  return result;
}
