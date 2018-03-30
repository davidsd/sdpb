//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "SDP_Solver_Parameters.hxx"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

int solve(const std::vector<boost::filesystem::path> &sdp_files,
          const boost::filesystem::path &out_file,
          const boost::filesystem::path &checkpoint_file_in,
          const boost::filesystem::path &checkpoint_file_out,
          SDP_Solver_Parameters parameters);

int main(int argc, char **argv)
{
  std::vector<boost::filesystem::path> sdp_files;
  boost::filesystem::path out_file;
  boost::filesystem::path checkpoint_file_in;
  boost::filesystem::path checkpoint_file_out;
  boost::filesystem::path param_file;

  SDP_Solver_Parameters parameters;

  po::options_description basic_options("Basic options");
  basic_options.add_options()("help,h", "Show this helpful message.")(
  "sdpFile,s",
  po::value<std::vector<boost::filesystem::path>>(&sdp_files)->required(),
  "SDP data file(s) in XML format. Use this option repeatedly to specify "
  "multiple data files.")(
  "paramFile,p", po::value<boost::filesystem::path>(&param_file),
  "Any parameter can optionally be set via this file in key=value "
  "format. Command line arguments override values in the parameter "
  "file.")("outFile,o", po::value<boost::filesystem::path>(&out_file),
           "The optimal solution is saved to this file in Mathematica "
           "format. Defaults to sdpFile with '.out' extension.")(
  "checkpointFile,c", po::value<boost::filesystem::path>(&checkpoint_file_out),
  "Checkpoints are saved to this file every checkpointInterval. Defaults "
  "to sdpFile with '.ck' extension.")(
  "initialCheckpointFile,i",
  po::value<boost::filesystem::path>(&checkpoint_file_in),
  "The initial checkpoint to load. Defaults to checkpointFile.");

  po::options_description solver_params_options("Solver parameters");
  solver_params_options.add_options()(
  "precision", po::value<int>(&parameters.precision)->default_value(400),
  "Precision in binary digits.  GMP will round up to the nearest "
  "multiple of 64 (or 32 on older systems).")(
  "maxThreads", po::value<int>(&parameters.maxThreads)->default_value(4),
  "Maximum number of threads to use for parallel calculation.")(
  "checkpointInterval",
  po::value<int>(&parameters.checkpointInterval)->default_value(3600),
  "Save checkpoints to checkpointFile every checkpointInterval seconds.")(
  "noFinalCheckpoint",
  po::bool_switch(&parameters.noFinalCheckpoint)->default_value(false),
  "Don't save a final checkpoint after terminating (useful when debugging).")(
  "findPrimalFeasible",
  po::bool_switch(&parameters.findPrimalFeasible)->default_value(false),
  "Terminate once a primal feasible solution is found.")(
  "findDualFeasible",
  po::bool_switch(&parameters.findDualFeasible)->default_value(false),
  "Terminate once a dual feasible solution is found.")(
  "detectPrimalFeasibleJump",
  po::bool_switch(&parameters.detectPrimalFeasibleJump)->default_value(false),
  "Terminate if a primal-step of 1 is taken. This often indicates that a "
  "primal feasible solution would be found if the precision were high "
  "enough. Try increasing either primalErrorThreshold or precision "
  "and run from the latest checkpoint.")(
  "detectDualFeasibleJump",
  po::bool_switch(&parameters.detectDualFeasibleJump)->default_value(false),
  "Terminate if a dual-step of 1 is taken. This often indicates that a "
  "dual feasible solution would be found if the precision were high "
  "enough. Try increasing either dualErrorThreshold or precision "
  "and run from the latest checkpoint.")(
  "maxIterations",
  po::value<int>(&parameters.maxIterations)->default_value(500),
  "Maximum number of iterations to run the solver.")(
  "maxRuntime", po::value<int>(&parameters.maxRuntime)->default_value(86400),
  "Maximum amount of time to run the solver in seconds.")(
  "dualityGapThreshold",
  po::value<Real>(&parameters.dualityGapThreshold)
  ->default_value(Real("1e-30")),
  "Threshold for duality gap (roughly the difference in primal and dual "
  "objective) at which the solution is considered "
  "optimal. Corresponds to SDPA's epsilonStar.")(
  "primalErrorThreshold",
  po::value<Real>(&parameters.primalErrorThreshold)
  ->default_value(Real("1e-30")),
  "Threshold for feasibility of the primal problem. Corresponds to SDPA's "
  "epsilonBar.")(
  "dualErrorThreshold",
  po::value<Real>(&parameters.dualErrorThreshold)->default_value(Real("1e-30")),
  "Threshold for feasibility of the dual problem. Corresponds to SDPA's "
  "epsilonBar.")(
  "initialMatrixScalePrimal",
  po::value<Real>(&parameters.initialMatrixScalePrimal)
  ->default_value(Real("1e20")),
  "The primal matrix X begins at initialMatrixScalePrimal times the "
  "identity matrix. Corresponds to SDPA's lambdaStar.")(
  "initialMatrixScaleDual",
  po::value<Real>(&parameters.initialMatrixScaleDual)
  ->default_value(Real("1e20")),
  "The dual matrix Y begins at initialMatrixScaleDual times the "
  "identity matrix. Corresponds to SDPA's lambdaStar.")(
  "feasibleCenteringParameter",
  po::value<Real>(&parameters.feasibleCenteringParameter)
  ->default_value(Real("0.1")),
  "Shrink the complementarity X Y by this factor when the primal and dual "
  "problems are feasible. Corresponds to SDPA's betaStar.")(
  "infeasibleCenteringParameter",
  po::value<Real>(&parameters.infeasibleCenteringParameter)
  ->default_value(Real("0.3")),
  "Shrink the complementarity X Y by this factor when either the primal "
  "or dual problems are infeasible. Corresponds to SDPA's betaBar.")(
  "stepLengthReduction",
  po::value<Real>(&parameters.stepLengthReduction)->default_value(Real("0.7")),
  "Shrink each newton step by this factor (smaller means slower, more "
  "stable convergence). Corresponds to SDPA's gammaStar.")(
  "choleskyStabilizeThreshold",
  po::value<Real>(&parameters.choleskyStabilizeThreshold)
  ->default_value(Real("1e-40")),
  "Adds stabilizing terms to the cholesky decomposition of the schur "
  "complement "
  "matrix for diagonal entries which are smaller than this threshold times "
  "the "
  "geometric mean of other diagonal entries. Somewhat higher "
  "choleskyStabilizeThreshold "
  "can improve numerical stability but if the threshold is large enough that "
  "a high "
  "proportion of eigenvalues are being stabilized, the computation will slow "
  "substantially.")(
  "maxComplementarity",
  po::value<Real>(&parameters.maxComplementarity)->default_value(Real("1e100")),
  "Terminate if the complementarity mu = Tr(X Y)/dim(X) exceeds this value.");

  po::options_description cmd_line_options;
  cmd_line_options.add(basic_options).add(solver_params_options);

  po::variables_map variables_map;

  try
    {
      po::store(po::parse_command_line(argc, argv, cmd_line_options),
                variables_map);

      if(variables_map.count("help"))
        {
          std::cout << cmd_line_options << '\n';
          return 0;
        }

      if(variables_map.count("paramFile"))
        {
          param_file = variables_map["paramFile"].as<boost::filesystem::path>();
          std::ifstream ifs(param_file.string().c_str());
          po::store(po::parse_config_file(ifs, solver_params_options),
                    variables_map);
        }

      po::notify(variables_map);

      if(!variables_map.count("outFile"))
        {
          out_file = sdp_files[0];
          out_file.replace_extension("out");
        }

      if(!variables_map.count("checkpointFile"))
        {
          checkpoint_file_out = sdp_files[0];
          checkpoint_file_out.replace_extension("ck");
        }

      if(!variables_map.count("initialCheckpointFile"))
        {
          checkpoint_file_in = checkpoint_file_out;
        }

      std::ofstream ofs(out_file.string().c_str());
      ofs.close();
      if(!ofs)
        {
          std::cerr << "Cannot write to outFile." << '\n';
          return 1;
        }
    }
  catch(po::error &e)
    {
      std::cerr << "ERROR: " << e.what() << '\n';
      std::cerr << cmd_line_options << '\n';
      return 1;
    }

  return solve(sdp_files, out_file, checkpoint_file_in, checkpoint_file_out,
               parameters);
}
