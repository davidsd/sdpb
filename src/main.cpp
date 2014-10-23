#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include "omp.h"
#include "boost/filesystem.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/program_options.hpp"
#include "types.h"
#include "Timers.h"
#include "SDP.h"
#include "parse.h"
#include "tests.h"
#include "SDPSolver.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setfill;
using std::setw;

using boost::filesystem::path;
using boost::posix_time::second_clock;

namespace po = boost::program_options;

Timers timers;

int solveSDP(const path &sdpFile,
             const path &outFile,
             const path &checkpointFile,
             SDPSolverParameters parameters) {

  setDefaultPrecision(parameters.precision);
  cout.precision(min(int(parameters.precision * 0.30102999566398114 + 5), 30));
  // Ensure all the Real parameters have the appropriate precision
  parameters.resetPrecision();
  omp_set_num_threads(parameters.maxThreads);

  cout << "SDPB started at " << second_clock::local_time() << endl;
  cout << "SDP file        : " << sdpFile        << endl;
  cout << "out file        : " << outFile        << endl;
  cout << "checkpoint file : " << checkpointFile << endl;
  cout << "using " << omp_get_max_threads() << " threads." << endl;

  cout << "\nParameters:\n";
  cout << parameters << endl;

  const SDP sdp = readBootstrapSDP(sdpFile);
  SDPSolver solver(sdp);

  if (exists(checkpointFile))
    solver.loadCheckpoint(checkpointFile);
  else
    solver.initialize(parameters);

  timers["Run solver"].start();
  timers["Save checkpoint"].start();
  printSolverHeader();
  SDPSolverTerminateReason reason = solver.run(parameters, checkpointFile);
  timers["Run solver"].stop();
  timers["Save checkpoint"].stop();

  cout << "-----" << setfill('-') << setw(100) << std::left << reason << endl << endl;
  cout << solver.status << endl;
  cout << timers << endl;

  if (!parameters.noFinalCheckpoint)
    solver.saveCheckpoint(checkpointFile);
  solver.saveSolution(reason, outFile);

  return 0;
}

int main(int argc, char** argv) {

  path sdpFile;
  path outFile;
  path checkpointFile;
  path paramFile;

  SDPSolverParameters parameters;

  po::options_description basicOptions("Basic options");
  basicOptions.add_options()
    ("help,h", "Show this helpful message.")
    ("sdpFile,s",
     po::value<path>(&sdpFile)->required(),
     "SDP data file in XML format.")
    ("paramFile,p",
     po::value<path>(&paramFile),
     "Any parameter can optionally be set via this file in key=value "
     "format. Command line arguments override values in the parameter "
     "file.")
    ("outFile,o",
     po::value<path>(&outFile),
     "The optimal solution is saved to this file in Mathematica "
     "format. Defaults to sdpFile with '.out' extension.")
    ("checkpointFile,c",
     po::value<path>(&checkpointFile),
     "Checkpoints are saved to this file every checkpointInterval. Defaults "
     "to sdpFile with '.ck' extension.")
    ;

  po::options_description solverParamsOptions("Solver parameters");
  solverParamsOptions.add_options()
    ("precision",
     po::value<int>(&parameters.precision)->default_value(400),
     "Precision in binary digits.  GMP will round up to the nearest "
     "multiple of 64.")
    ("maxThreads",
     po::value<int>(&parameters.maxThreads)->default_value(4),
     "Maximum number of threads to use for parallel calculation.")
    ("checkpointInterval",
     po::value<int>(&parameters.checkpointInterval)->default_value(3600),
     "Save checkpoints to checkpointFile every checkpointInterval seconds.")
    ("noFinalCheckpoint",
     po::bool_switch(&parameters.noFinalCheckpoint)->default_value(false),
     "Don't save a final checkpoint after terminating (useful when debugging).")
    ("findDualFeasible",
     po::bool_switch(&parameters.findDualFeasible)->default_value(false),
     "Terminate once a dual feasible solution is found.")
    ("detectDualFeasibleJump",
     po::bool_switch(&parameters.detectDualFeasibleJump)->default_value(false),
     "Terminate if a dual-step of 1 is taken. This often indicates that a "
     "dual feasible solution would be found if the precision were high "
     "enough. Try increasing either dualErrorThreshold or precision "
     "and run from the latest checkpoint.")
    ("maxIterations",
     po::value<int>(&parameters.maxIterations)->default_value(500),
     "Maximum number of iterations to run the solver.")
    ("maxRuntime",
     po::value<int>(&parameters.maxRuntime)->default_value(86400),
     "Maximum amount of time to run the solver in seconds.")
    ("dualityGapThreshold",
     po::value<Real>(&parameters.dualityGapThreshold)->default_value(Real("1e-30")),
     "Threshold for duality gap (roughly the difference in primal and dual "
     "objective) at which the solution is considered "
     "optimal. Corresponds to SDPA's epsilonStar.")
    ("primalErrorThreshold",
     po::value<Real>(&parameters.primalErrorThreshold)->default_value(Real("1e-30")),
     "Threshold for feasibility of the primal problem. Corresponds to SDPA's "
     "epsilonBar.")
    ("dualErrorThreshold",
     po::value<Real>(&parameters.dualErrorThreshold)->default_value(Real("1e-30")),
     "Threshold for feasibility of the dual problem. Corresponds to SDPA's epsilonBar.")
    ("initialMatrixScalePrimal",
     po::value<Real>(&parameters.initialMatrixScalePrimal)->default_value(Real("1e20")),
     "The primal matrix X begins at initialMatrixScalePrimal times the "
     "identity matrix. Corresponds to SDPA's lambdaStar.")
    ("initialMatrixScaleDual",
     po::value<Real>(&parameters.initialMatrixScaleDual)->default_value(Real("1e20")),
     "The dual matrix Y begins at initialMatrixScaleDual times the "
     "identity matrix. Corresponds to SDPA's lambdaStar.")
    ("feasibleCenteringParameter",
     po::value<Real>(&parameters.feasibleCenteringParameter)->default_value(Real("0.1")),
     "Shrink the complementarity X Y by this factor when the primal and dual "
     "problems are feasible. Corresponds to SDPA's betaStar.")
    ("infeasibleCenteringParameter",
     po::value<Real>(&parameters.infeasibleCenteringParameter)->default_value(Real("0.3")),
     "Shrink the complementarity X Y by this factor when either the primal "
     "or dual problems are infeasible. Corresponds to SDPA's betaBar.")
    ("stepLengthReduction",
     po::value<Real>(&parameters.stepLengthReduction)->default_value(Real("0.7")),
     "Shrink each newton step by this factor (smaller means slower, more "
     "stable convergence). Corresponds to SDPA's gammaStar.")
    ("maxComplementarity",
     po::value<Real>(&parameters.maxComplementarity)->default_value(Real("1e100")),
     "Terminate if the complementarity mu = Tr(X Y)/dim(X) exceeds this value.")
    ;
    
  po::options_description cmdLineOptions;
  cmdLineOptions.add(basicOptions).add(solverParamsOptions);

  po::variables_map variablesMap;

  try {
    po::store(po::parse_command_line(argc, argv, cmdLineOptions), variablesMap);

    if (variablesMap.count("help")) {
      cout << cmdLineOptions << endl;
      return 0;
    }

    if (variablesMap.count("paramFile")) {
      paramFile = variablesMap["paramFile"].as<path>();
      std::ifstream ifs(paramFile.string().c_str());
      po::store(po::parse_config_file(ifs, solverParamsOptions), variablesMap);
    }

    po::notify(variablesMap);

    if (!variablesMap.count("outFile")) {
      outFile = sdpFile;
      outFile.replace_extension("out");
    }

    if (!variablesMap.count("checkpointFile")) {
      checkpointFile = sdpFile;
      checkpointFile.replace_extension("ck");
    }
  } catch(po::error& e) {
    cerr << "ERROR: " << e.what() << endl;
    cerr << cmdLineOptions << endl; 
    return 1; 
  } 

  return solveSDP(sdpFile, outFile, checkpointFile, parameters);
}
