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
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include "../types.hxx"
#include "../Timers.hxx"
#include "../SDP.hxx"
#include "../read_bootstrap_sdp.hxx"
#include "../SDPSolver.hxx"

using std::cout;
using std::cerr;
using std::endl;
using std::setfill;
using std::setw;

using boost::filesystem::path;
using boost::posix_time::second_clock;

namespace po = boost::program_options;

Timers timers;

int solve(const std::vector<boost::filesystem::path> &sdpFiles,
          const boost::filesystem::path &outFile,
          const boost::filesystem::path &checkpointFileIn,
          const boost::filesystem::path &checkpointFileOut,
          SDPSolverParameters parameters) {
  // Set the default precision of all Real numbers to that specified
  // by the 'precision' parameter.
  mpf_set_default_prec(parameters.precision);

  // Set cout to print the appropriate number of digits
  cout.precision(min(static_cast<int>(parameters.precision * 0.31 + 5), 30));

  // Ensure all the Real parameters have the appropriate precision
  parameters.resetPrecision();


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
  SDPSolver solver(read_bootstrap_sdp(sdpFiles), parameters);
  
  if (exists(checkpointFileIn))
    solver.loadCheckpoint(checkpointFileIn);

  timers["Solver runtime"].start();
  timers["Last checkpoint"].start();
  SDPSolverTerminateReason reason = solver.run(checkpointFileOut);
  timers["Solver runtime"].stop();
  //SDPSolverTerminateReason reason;
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

  timers.writeMFile( outFile.string()+".profiling" );

  return 0;
}
