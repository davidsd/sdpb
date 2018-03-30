//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "SDP.hxx"
#include "SDP_Solver.hxx"
#include "../Timers.hxx"

#include <algorithm>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>

// FIXME: Pass this around instead of having a global.
Timers timers;

SDP read_bootstrap_sdp(const std::vector<boost::filesystem::path> &sdp_files);


int solve(const std::vector<boost::filesystem::path> &sdpFiles,
          const boost::filesystem::path &outFile,
          const boost::filesystem::path &checkpointFileIn,
          const boost::filesystem::path &checkpointFileOut,
          SDP_Solver_Parameters parameters)
{
  // Set the default precision of all Real numbers to that specified
  // by the 'precision' parameter.
  mpf_set_default_prec(parameters.precision);

  // Set std::cout to print the appropriate number of digits
  std::cout.precision(
  min(static_cast<int>(parameters.precision * 0.31 + 5), 30));

  // Ensure all the Real parameters have the appropriate precision
  parameters.resetPrecision();

  std::cout << "SDPB started at "
            << boost::posix_time::second_clock::local_time() << '\n';
  for(auto const &sdpFile : sdpFiles)
    { std::cout << "SDP file        : " << sdpFile << '\n'; }
    std::cout << "out file        : " << outFile << '\n';
  std::cout << "checkpoint in   : " << checkpointFileIn << '\n';
  std::cout << "checkpoint out  : " << checkpointFileOut << '\n';

  std::cout << "\nParameters:\n";
  std::cout << parameters << '\n';

  // Read an SDP from sdpFile and create a solver for it
  SDP_Solver solver(read_bootstrap_sdp(sdpFiles), parameters);

  if(exists(checkpointFileIn))
    solver.loadCheckpoint(checkpointFileIn);

  timers["Solver runtime"].start();
  timers["Last checkpoint"].start();
  SDPSolverTerminateReason reason = solver.run(checkpointFileOut);
  timers["Solver runtime"].stop();
  std::cout << "-----" << std::setfill('-') << std::setw(116) << std::left
            << reason << '\n';
  std::cout << '\n';
  std::cout << "primalObjective = " << solver.primalObjective << '\n';
  std::cout << "dualObjective   = " << solver.dualObjective << '\n';
  std::cout << "dualityGap      = " << solver.dualityGap << '\n';
  std::cout << "primalError     = " << solver.primalError << '\n';
  std::cout << "dualError       = " << solver.dualError << '\n';
  std::cout << '\n';

  if(!parameters.noFinalCheckpoint)
    solver.saveCheckpoint(checkpointFileOut);
  timers["Last checkpoint"].stop();
  solver.saveSolution(reason, outFile);

  std::cout << '\n' << timers;

  timers.writeMFile(outFile.string() + ".profiling");

  return 0;
}
