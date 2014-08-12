#include <iostream>
#include <ostream>
#include "omp.h"
#include "boost/archive/text_oarchive.hpp"
#include "boost/archive/text_iarchive.hpp"
#include "boost/filesystem.hpp"
#include "boost/filesystem/fstream.hpp"
#include "SDPSolver.h"
#include "serialize.h"
#include "Timers.h"

using boost::filesystem::path;
using std::cout;

void printSolverHeader() {
  cout << "\n     time    mu        P-obj       D-obj      gap         P-err         D-err       P-step   D-step   beta\n";
  cout << "-----------------------------------------------------------------------------------------------------------\n";
}

void printSolverInfo(int iteration,
                     Real mu,
                     SDPSolverStatus status,
                     bool isPrimalFeasible,
                     bool isDualFeasible,
                     Real primalStepLength,
                     Real dualStepLength,
                     Real betaCorrector) {
  Real time = Real(timers["Run solver"].elapsed().wall)/1000000000;
  gmp_fprintf(stdout,
              "%3d  %-7.3Fg %-8.1Fe %-+11.2Fe %-+11.2Fe %-9.2Fe %s%-+10.2Fe%s  %s%-+10.2Fe%s  %-8.3Fg %-8.3Fg %-4.2Fg\n",
              iteration,
              time.get_mpf_t(),
              mu.get_mpf_t(),
              status.primalObjective.get_mpf_t(),
              status.dualObjective.get_mpf_t(),
              status.dualityGap().get_mpf_t(),
              isPrimalFeasible ? "|" : " ", status.primalError.get_mpf_t(), isPrimalFeasible ? "|" : " ",
              isDualFeasible   ? "|" : " ", status.dualError.get_mpf_t(),   isDualFeasible   ? "|" : " ",
              primalStepLength.get_mpf_t(),
              dualStepLength.get_mpf_t(),
              betaCorrector.get_mpf_t());
}

ostream& operator<<(ostream& os, const SDPSolverParameters& p) {
  os << "maxIterations                = " << p.maxIterations                << endl;
  os << "maxRuntime                   = " << p.maxRuntime                   << endl;
  os << "checkpointInterval           = " << p.checkpointInterval           << endl;
  os << "precision(actual)            = " << p.precision << "(" << mpf_get_default_prec() << ")" << endl;
  os << "maxThreads                   = " << p.maxThreads                   << endl;
  os << "dualityGapThreshold          = " << p.dualityGapThreshold          << endl;
  os << "primalErrorThreshold         = " << p.primalErrorThreshold         << endl;
  os << "dualErrorThreshold           = " << p.dualErrorThreshold           << endl;
  os << "initialMatrixScale           = " << p.initialMatrixScale           << endl;
  os << "feasibleCenteringParameter   = " << p.feasibleCenteringParameter   << endl;
  os << "infeasibleCenteringParameter = " << p.infeasibleCenteringParameter << endl;
  os << "stepLengthReduction          = " << p.stepLengthReduction          << endl;
  os << "maxDualObjective             = " << p.maxDualObjective             << endl;
  return os;
}

ostream &operator<<(ostream& os, const SDPSolverTerminateReason& r) {
  switch(r) {
  case PrimalDualOptimal:
    os << "found primal-dual optimal solution.";
    break;
  case MaxIterationsExceeded:
    os << "maxIterations exceeded.";
    break;
  case MaxRuntimeExceeded:
    os << "maxRuntime exceeded.";
    break;
  case DualFeasibleMaxObjectiveExceeded:
    os << "found dual feasible solution with dualObjective exceeding maxDualObjective.";
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
  boost::serialization::serializeSDPSolverState(ar, x, X, Y);
}

void SDPSolver::loadCheckpoint(const path &checkpointFile) {
  boost::filesystem::ifstream ifs(checkpointFile);
  boost::archive::text_iarchive ar(ifs);
  cout << "Loading checkpoint from : " << checkpointFile << endl;
  boost::serialization::serializeSDPSolverState(ar, x, X, Y);
}

void SDPSolver::initialize(const SDPSolverParameters &parameters) {
  fillVector(x, 0);
  X.setZero();
  X.addDiagonal(parameters.initialMatrixScale);
  Y.setZero();
  Y.addDiagonal(parameters.initialMatrixScale);
}

void SDPSolver::saveSolution(const path &outFile) {
  boost::filesystem::ofstream ofs(outFile);
  cout << "Saving solution to: " << outFile << endl;
  ofs.precision(int(status.primalObjective.get_prec() * 0.30102999566398114 + 5));
  ofs << "primalObjective = " << status.primalObjective << ";\n";
  ofs << "dualObjective   = " << status.dualObjective   << ";\n";
  ofs << "dualityGap      = " << status.dualityGap()    << ";\n";
  ofs << "primalError     = " << status.primalError     << ";\n";
  ofs << "dualError       = " << status.dualError       << ";\n";
  ofs << "runtime         = " << float(timers["Run solver"].elapsed().wall)/1000000000 << ";\n";
  ofs << "freeVariables   = " << freeVariableSolution() << ";\n";
  ofs << "x = " << x << ";\n";
  ofs << "X = " << X << ";\n";
  ofs << "Y = " << Y << ";\n";
}
