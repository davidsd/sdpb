//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================


#include <iostream>
#include <ostream>
#include <sstream>
#include "omp.h"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
//Tweak to allow Ubuntu-14.04/gcc-4.8.4 and similar environments to compile
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS
#include "boost/filesystem/fstream.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "SDPSolver.h"
#include "serialize.h"
#include "Timers.h"

using boost::filesystem::path;
using boost::posix_time::time_duration;
using boost::posix_time::microseconds;
using std::cout;

ostream& operator<<(ostream& os, const SDPSolverParameters& p) {
  os << std::boolalpha;
  os << "maxIterations                = " << p.maxIterations                << endl;
  os << "maxRuntime                   = " << p.maxRuntime                   << endl;
  os << "checkpointInterval           = " << p.checkpointInterval           << endl;
  os << "noFinalCheckpoint            = " << p.noFinalCheckpoint            << endl;
  os << "findPrimalFeasible           = " << p.findPrimalFeasible           << endl;
  os << "findDualFeasible             = " << p.findDualFeasible             << endl;
  os << "detectPrimalFeasibleJump     = " << p.detectPrimalFeasibleJump     << endl;
  os << "detectDualFeasibleJump       = " << p.detectDualFeasibleJump       << endl;
  os << "precision(actual)            = " << p.precision << "(" << mpf_get_default_prec() << ")" << endl;
  os << "maxThreads(using)            = " << p.maxThreads << "(" << omp_get_max_threads() << ")" << endl;
  os << "dualityGapThreshold          = " << p.dualityGapThreshold          << endl;
  os << "primalErrorThreshold         = " << p.primalErrorThreshold         << endl;
  os << "dualErrorThreshold           = " << p.dualErrorThreshold           << endl;
  os << "initialMatrixScalePrimal     = " << p.initialMatrixScalePrimal     << endl;
  os << "initialMatrixScaleDual       = " << p.initialMatrixScaleDual       << endl;
  os << "feasibleCenteringParameter   = " << p.feasibleCenteringParameter   << endl;
  os << "infeasibleCenteringParameter = " << p.infeasibleCenteringParameter << endl;
  os << "stepLengthReduction          = " << p.stepLengthReduction          << endl;
  os << "choleskyStabilizeThreshold   = " << p.choleskyStabilizeThreshold   << endl;
  os << "maxComplementarity           = " << p.maxComplementarity           << endl;
  return os;
}

ostream &operator<<(ostream& os, const SDPSolverTerminateReason& r) {
  switch (r) {
  case PrimalDualOptimal:
    os << "found primal-dual optimal solution";
    break;
  case PrimalFeasible:
    os << "found primal feasible solution";
    break;
  case DualFeasible:
    os << "found dual feasible solution";
    break;
  case PrimalFeasibleJumpDetected:
    os << "primal feasible jump detected";
    break;
  case DualFeasibleJumpDetected:
    os << "dual feasible jump detected";
    break;
  case MaxIterationsExceeded:
    os << "maxIterations exceeded";
    break;
  case MaxRuntimeExceeded:
    os << "maxRuntime exceeded";
    break;
  case MaxComplementarityExceeded:
    os << "maxComplementarity exceeded";
    break;
  }
  return os;
}

void SDPSolver::printHeader() {
  cout << "\n     time      mu        P-obj       D-obj      gap         P-err       D-err      P-step   D-step   beta  dim/stabilized\n";
  cout << "-------------------------------------------------------------------------------------------------------------------------\n";
}

void SDPSolver::printIteration(int iteration,
                               Real mu,
                               Real primalStepLength,
                               Real dualStepLength,
                               Real betaCorrector) {
  time_duration td(microseconds(timers["Solver runtime"].elapsed().wall)/1000);
  std::stringstream ss;
  ss << td;
  gmp_fprintf(stdout,
              "%3d  %s  %-8.1Fe %-+11.2Fe %-+11.2Fe %-9.2Fe  %-+10.2Fe  %-+10.2Fe  %-8.3Fg %-8.3Fg %-4.2Fg  %d/%d",
              iteration,
              ss.str().substr(0, 8).c_str(),
              mu.get_mpf_t(),
              primalObjective.get_mpf_t(),
              dualObjective.get_mpf_t(),
              dualityGap.get_mpf_t(),
              primalError.get_mpf_t(),
              dualError.get_mpf_t(),
              primalStepLength.get_mpf_t(),
              dualStepLength.get_mpf_t(),
              betaCorrector.get_mpf_t(),
              static_cast<int>(sdp.dualObjective.size()),
              Q.rows);
  cout << endl;
}

void backupCheckpointFile(path const& checkpointFile) {
  path backupFile(checkpointFile);
  backupFile.replace_extension(".ck.bk");
  cout << "Backing up checkpoint to: " << backupFile << endl;
  copy_file(checkpointFile, backupFile, boost::filesystem::copy_option::overwrite_if_exists);
}

void SDPSolver::saveCheckpoint(const path &checkpointFile) {
  if (exists(checkpointFile))
    backupCheckpointFile(checkpointFile);
  boost::filesystem::ofstream ofs(checkpointFile);
  boost::archive::text_oarchive ar(ofs);
  cout << "Saving checkpoint to    : " << checkpointFile << endl;
  boost::serialization::serializeSDPSolverState(ar, x, X, y, Y);
  timers["Last checkpoint"].start();
}

void SDPSolver::loadCheckpoint(const path &checkpointFile) {
  boost::filesystem::ifstream ifs(checkpointFile);
  boost::archive::text_iarchive ar(ifs);
  cout << "Loading checkpoint from : " << checkpointFile << endl;
  boost::serialization::serializeSDPSolverState(ar, x, X, y, Y);
}

void SDPSolver::saveSolution(const SDPSolverTerminateReason terminateReason, const path &outFile) {
  boost::filesystem::ofstream ofs(outFile);
  float runtime = static_cast<float>(timers["Solver runtime"].elapsed().wall)/1000000000;
  cout << "Saving solution to      : " << outFile << endl;
  ofs.precision(static_cast<int>(primalObjective.get_prec() * 0.31 + 5));
  ofs << "terminateReason = \"" << terminateReason << "\";\n";
  ofs << "primalObjective = " << primalObjective   << ";\n";
  ofs << "dualObjective   = " << dualObjective     << ";\n";
  ofs << "dualityGap      = " << dualityGap        << ";\n";
  ofs << "primalError     = " << primalError       << ";\n";
  ofs << "dualError       = " << dualError         << ";\n";
  ofs << "runtime         = " << runtime           << ";\n";
  ofs << "y = " << y << ";\n";
  // ofs << "Y = " << Y << ";\n";
  ofs << "x = " << x << ";\n";
  // ofs << "X = " << X << ";\n";
}
