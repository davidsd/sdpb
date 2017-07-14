//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include "omp.h"
//Tweak to allow Ubuntu-14.04/gcc-4.8.4 and similar environments to compile
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include "boost/filesystem.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/program_options.hpp"
#include "types.h"
#include "Timers.h"
#include "SDP.h"
#include "parse.h"
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

int solveSDP(const vector<path> &sdpFiles,
             const path &outFile,
             const path &checkpointFileIn,
             const path &checkpointFileOut,
             SDPSolverParameters parameters) {
  // Set the default precision of all Real numbers to that specified
  // by the 'precision' parameter.
  mpf_set_default_prec(parameters.precision);

  // Set cout to print the appropriate number of digits
  cout.precision(min(static_cast<int>(parameters.precision * 0.31 + 5), 30));

  // Ensure all the Real parameters have the appropriate precision
  parameters.resetPrecision();

  // Set the maximum number of threads for parallel loops
  omp_set_num_threads(parameters.maxThreads);

  cout << "SDPB started at " << second_clock::local_time() << endl;
  for (auto const& sdpFile: sdpFiles) {
    cout << "SDP file        : " << sdpFile        << endl;
  }
  cout << "out file        : " << outFile        << endl;
  cout << "checkpoint in   : " << checkpointFileIn << endl;
  cout << "checkpoint out  : " << checkpointFileOut << endl;

  cout << "\nParameters:\n";
  cout << parameters << endl;

  // Read an SDP from sdpFile and create a solver for it
  SDPSolver solver(readBootstrapSDP(sdpFiles), parameters);

  if (exists(checkpointFileIn))
    solver.loadCheckpoint(checkpointFileIn);

  timers["Solver runtime"].start();
  timers["Last checkpoint"].start();
  SDPSolverTerminateReason reason = solver.run(checkpointFileOut);
  timers["Solver runtime"].stop();

  cout << "-----" << setfill('-') << setw(116) << std::left << reason << endl;
  cout << endl;
  cout << "primalObjective = " << solver.primalObjective << endl;
  cout << "dualObjective   = " << solver.dualObjective   << endl;
  cout << "dualityGap      = " << solver.dualityGap      << endl;
  cout << "primalError     = " << solver.primalError     << endl;
  cout << "dualError       = " << solver.dualError       << endl;
  cout << endl;

  if (!parameters.noFinalCheckpoint)
    solver.saveCheckpoint(checkpointFileOut);
  timers["Last checkpoint"].stop();
  solver.saveSolution(reason, outFile);

  cout << endl << timers;

  return 0;
}

int main(int argc, char** argv) {
  vector<path> sdpFiles;
  path outFile;
  path checkpointFileIn;
  path checkpointFileOut;
  path paramFile;

  SDPSolverParameters parameters;

  po::options_description basicOptions("Basic options");
  basicOptions.add_options()
    ("help,h", "Show this helpful message.")
    ("sdpFile,s",
     po::value< vector<path> >(&sdpFiles)->required(),
     "SDP data file(s) in XML format. Use this option repeatedly to specify multiple data files.")
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
     po::value<path>(&checkpointFileOut),
     "Checkpoints are saved to this file every checkpointInterval. Defaults "
     "to sdpFile with '.ck' extension.")
    ("initialCheckpointFile,i",
     po::value<path>(&checkpointFileIn),
     "The initial checkpoint to load. Defaults to checkpointFile.")
    ;

  po::options_description solverParamsOptions("Solver parameters");
  solverParamsOptions.add_options()
    ("precision",
     po::value<int>(&parameters.precision)->default_value(400),
     "Precision in binary digits.  GMP will round up to the nearest "
     "multiple of 64 (or 32 on older systems).")
    ("maxThreads",
     po::value<int>(&parameters.maxThreads)->default_value(4),
     "Maximum number of threads to use for parallel calculation.")
    ("checkpointInterval",
     po::value<int>(&parameters.checkpointInterval)->default_value(3600),
     "Save checkpoints to checkpointFile every checkpointInterval seconds.")
    ("noFinalCheckpoint",
     po::bool_switch(&parameters.noFinalCheckpoint)->default_value(false),
     "Don't save a final checkpoint after terminating (useful when debugging).")
    ("findPrimalFeasible",
     po::bool_switch(&parameters.findPrimalFeasible)->default_value(false),
     "Terminate once a primal feasible solution is found.")
    ("findDualFeasible",
     po::bool_switch(&parameters.findDualFeasible)->default_value(false),
     "Terminate once a dual feasible solution is found.")
    ("detectPrimalFeasibleJump",
     po::bool_switch(&parameters.detectPrimalFeasibleJump)->default_value(false),
     "Terminate if a primal-step of 1 is taken. This often indicates that a "
     "primal feasible solution would be found if the precision were high "
     "enough. Try increasing either primalErrorThreshold or precision "
     "and run from the latest checkpoint.")
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
    ("choleskyStabilizeThreshold",
     po::value<Real>(&parameters.choleskyStabilizeThreshold)->default_value(Real("1e-40")),
     "Adds stabilizing terms to the cholesky decomposition of the schur complement "
     "matrix for diagonal entries which are smaller than this threshold times the "
     "geometric mean of other diagonal entries. Somewhat higher choleskyStabilizeThreshold "
     "can improve numerical stability but if the threshold is large enough that a high "
     "proportion of eigenvalues are being stabilized, the computation will slow substantially.")
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
      outFile = sdpFiles[0];
      outFile.replace_extension("out");
    }

    if (!variablesMap.count("checkpointFile")) {
      checkpointFileOut = sdpFiles[0];
      checkpointFileOut.replace_extension("ck");
    }

    if (!variablesMap.count("initialCheckpointFile")) {
      checkpointFileIn = checkpointFileOut;
    }

    std::ofstream ofs(outFile.string().c_str());
    ofs.close();
    if (!ofs) {
      cerr << "Cannot write to outFile." << endl;
      return 1;
    }
  } catch(po::error& e) {
    cerr << "ERROR: " << e.what() << endl;
    cerr << cmdLineOptions << endl;
    return 1;
  }

  return solveSDP(sdpFiles, outFile, checkpointFileIn, checkpointFileOut, parameters);
}
