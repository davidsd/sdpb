#include <iostream>
#include <ostream>
#include <sstream>
#include "omp.h"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "SDPSolver.h"
#include "serialize.h"
#include "Timers.h"

using boost::filesystem::path;
using boost::posix_time::time_duration;
using boost::posix_time::microseconds;
using std::cout;

void printSolverHeader() {
  cout << "\n     time      mu        P-obj       D-obj      gap         P-err       D-err      P-step   D-step   beta  dim/stabilized\n";
  cout << "-------------------------------------------------------------------------------------------------------------------------\n";
}

void printSolverInfo(int iteration,
                     Real mu,
                     SDPSolverStatus status,
                     Real primalStepLength,
                     Real dualStepLength,
                     Real betaCorrector,
                     int dualObjectiveSize,
                     int Qrows) {
  Real time = Real(timers["Solver runtime"].elapsed().wall)/1000000000;
  time_duration td(microseconds(timers["Solver runtime"].elapsed().wall)/1000);
  std::stringstream ss;
  ss << td;
  gmp_fprintf(stdout,
              "%3d  %s  %-8.1Fe %-+11.2Fe %-+11.2Fe %-9.2Fe  %-+10.2Fe  %-+10.2Fe  %-8.3Fg %-8.3Fg %-4.2Fg  %d/%d",
              iteration,
              ss.str().substr(0,8).c_str(),
              mu.get_mpf_t(),
              status.primalObjective.get_mpf_t(),
              status.dualObjective.get_mpf_t(),
              status.dualityGap().get_mpf_t(),
              status.primalError.get_mpf_t(),
              status.dualError.get_mpf_t(),
              primalStepLength.get_mpf_t(),
              dualStepLength.get_mpf_t(),
              betaCorrector.get_mpf_t(),
              dualObjectiveSize,
              Qrows);
  cout << endl;
}

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
  switch(r) {
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

ostream& operator<<(ostream& os, const SDPSolverStatus& s) {
  os << "primalObjective = " << s.primalObjective << endl;
  os << "dualObjective   = " << s.dualObjective << endl;
  os << "dualityGap      = " << s.dualityGap() << endl;
  os << "primalError     = " << s.primalError << endl;
  os << "dualError       = " << s.dualError << endl;
  return os;
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
  cout << "Saving solution to      : " << outFile << endl;
  ofs.precision(int(status.primalObjective.get_prec() * 0.30102999566398114 + 5));
  ofs << "terminateReason = \"" << terminateReason      << "\";\n";
  ofs << "primalObjective = " << status.primalObjective << ";\n";
  ofs << "dualObjective   = " << status.dualObjective   << ";\n";
  ofs << "dualityGap      = " << status.dualityGap()    << ";\n";
  ofs << "primalError     = " << status.primalError     << ";\n";
  ofs << "dualError       = " << status.dualError       << ";\n";
  ofs << "runtime         = " << float(timers["Solver runtime"].elapsed().wall)/1000000000 << ";\n";
  ofs << "y = " << y << ";\n";
  ofs << "Y = " << Y << ";\n";
  ofs << "x = " << x << ";\n";
  ofs << "X = " << X << ";\n";
}
