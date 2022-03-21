#include "../Dynamical_Solver_Parameters.hxx"

#include <boost/program_options.hpp>
#include <vector>

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
   result.add_options()("prevTotalIterations", 
                              boost::program_options::value<size_t>(&total_iterations)->default_value(0),
                              "The number of total iterations finished before this solver starts to run. "
                              "The default value is set to 0.");
  result.add_options()("updateSdpThresholdMax",
                              boost::program_options::value<El::BigFloat>(&update_sdp_threshold_max)->default_value(1),
                              "Take a step in the external parameters, "
                              "that is to regenerate the sdp files if the step size is smaller than the threshold. "
                              "The default value is set to 1.");  
  result.add_options()("updateSdpThresholdMin",
                              boost::program_options::value<El::BigFloat>(&update_sdp_threshold_min)->default_value(0),
                              "Take a step in the external parameters, "
                              "that is to regenerate the sdp files if the step size is bigger than the threshold. "
                              "The default value is set to 0.");  
  result.add_options()("findBoundary", 
                              boost::program_options::value<bool>(&find_boundary) -> default_value(0), 
                              "True if the program is required by the user to find a point on the boundary of the island. " 
                              "The default value is set to False.");
  result.add_options()("findBoundaryObjThreshold",
                              boost::program_options::value<El::BigFloat>(&find_boundary_obj_threshold)->default_value(0),
                              "Continue to move towards the boundary if the primal and dual objectives are not sufficiently close to zero. " 
                              "The default value is set to 0.");
  result.add_options()("externalCoor", 
                              boost::program_options::value<std::vector<El::BigFloat>>(&external_coor)->multitoken(), 
                              "The values of the external variables that were used to produce the central sdp file.");
  result.add_options()("searchDirection", 
                              boost::program_options::value<std::vector<El::BigFloat>>(&search_direction)->multitoken(), 
                              "User-specified directional vector in which the program will looks for a zero. ");
  result.add_options()("lagMultiplier", 
                              boost::program_options::value<El::BigFloat>(&lag_multiplier_lambda), 
                              "The Lagrange multiplier used to find the boundary of an island. ");
  result.add_options()("muDirectionMode", 
                              boost::program_options::value<int>(&mu_last_direction),
                              "To decrease or increase mu based on the former step. ");

  result.add(solver_parameters.options()); 

  return result;
}
