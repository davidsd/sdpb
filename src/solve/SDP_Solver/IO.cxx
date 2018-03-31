//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "../SDP_Solver.hxx"
#include "../Timers.hxx"
#include "serialize.hxx"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <ostream>
#include <sstream>

std::ostream &operator<<(std::ostream &os, const SDP_Solver_Parameters &p)
{
  os << std::boolalpha;
  os << "maxIterations                = " << p.maxIterations << '\n';
  os << "maxRuntime                   = " << p.maxRuntime << '\n';
  os << "checkpointInterval           = " << p.checkpointInterval << '\n';
  os << "noFinalCheckpoint            = " << p.noFinalCheckpoint << '\n';
  os << "findPrimalFeasible           = " << p.findPrimalFeasible << '\n';
  os << "findDualFeasible             = " << p.findDualFeasible << '\n';
  os << "detectPrimalFeasibleJump     = " << p.detectPrimalFeasibleJump
     << '\n';
  os << "detectDualFeasibleJump       = " << p.detectDualFeasibleJump << '\n';
  os << "precision(actual)            = " << p.precision << "("
     << mpf_get_default_prec() << ")" << '\n';
  os << "dualityGapThreshold          = " << p.dualityGapThreshold << '\n';
  os << "primalErrorThreshold         = " << p.primalErrorThreshold << '\n';
  os << "dualErrorThreshold           = " << p.dualErrorThreshold << '\n';
  os << "initialMatrixScalePrimal     = " << p.initialMatrixScalePrimal
     << '\n';
  os << "initialMatrixScaleDual       = " << p.initialMatrixScaleDual << '\n';
  os << "feasibleCenteringParameter   = " << p.feasibleCenteringParameter
     << '\n';
  os << "infeasibleCenteringParameter = " << p.infeasibleCenteringParameter
     << '\n';
  os << "stepLengthReduction          = " << p.stepLengthReduction << '\n';
  os << "choleskyStabilizeThreshold   = " << p.choleskyStabilizeThreshold
     << '\n';
  os << "maxComplementarity           = " << p.maxComplementarity << '\n';
  return os;
}

std::ostream &
operator<<(std::ostream &os, const SDP_Solver_Terminate_Reason &r)
{
  switch(r)
    {
    case SDP_Solver_Terminate_Reason::PrimalDualOptimal:
      os << "found primal-dual optimal solution";
      break;
    case SDP_Solver_Terminate_Reason::PrimalFeasible:
      os << "found primal feasible solution";
      break;
    case SDP_Solver_Terminate_Reason::DualFeasible:
      os << "found dual feasible solution";
      break;
    case SDP_Solver_Terminate_Reason::PrimalFeasibleJumpDetected:
      os << "primal feasible jump detected";
      break;
    case SDP_Solver_Terminate_Reason::DualFeasibleJumpDetected:
      os << "dual feasible jump detected";
      break;
    case SDP_Solver_Terminate_Reason::MaxIterationsExceeded:
      os << "maxIterations exceeded";
      break;
    case SDP_Solver_Terminate_Reason::MaxRuntimeExceeded:
      os << "maxRuntime exceeded";
      break;
    case SDP_Solver_Terminate_Reason::MaxComplementarityExceeded:
      os << "maxComplementarity exceeded";
      break;
    }
  return os;
}

void SDP_Solver::printHeader()
{
  std::cout
    << "\n     time      mu        P-obj       D-obj      gap         "
       "P-err       D-err      P-step   D-step   beta  dim/stabilized\n";
  std::cout
    << "-------------------------------------------------------------------"
       "------------------------------------------------------\n";
}

void SDP_Solver::printIteration(int iteration, Real mu, Real primalStepLength,
                                Real dualStepLength, Real betaCorrector)
{
  boost::posix_time::time_duration td(
    boost::posix_time::microseconds(timers["Solver runtime"].elapsed().wall)
    / 1000);
  std::stringstream ss;
  ss << td;
  gmp_fprintf(stdout,
              "%3d  %s  %-8.1Fe %-+11.2Fe %-+11.2Fe %-9.2Fe  %-+10.2Fe  "
              "%-+10.2Fe  %-8.3Fg %-8.3Fg %-4.2Fg  %d/%d",
              iteration, ss.str().substr(0, 8).c_str(), mu.get_mpf_t(),
              primalObjective.get_mpf_t(), dualObjective.get_mpf_t(),
              dualityGap.get_mpf_t(), primalError.get_mpf_t(),
              dualError.get_mpf_t(), primalStepLength.get_mpf_t(),
              dualStepLength.get_mpf_t(), betaCorrector.get_mpf_t(),
              static_cast<int>(sdp.dualObjective.size()), Q.rows);
  std::cout << '\n';
}

void backupCheckpointFile(const boost::filesystem::path &checkpointFile)
{
  boost::filesystem::path backupFile(checkpointFile);
  backupFile.replace_extension(".ck.bk");
  std::cout << "Backing up checkpoint to: " << backupFile << '\n';
  copy_file(checkpointFile, backupFile,
            boost::filesystem::copy_option::overwrite_if_exists);
}

void SDP_Solver::saveCheckpoint(const boost::filesystem::path &checkpointFile)
{
  if(exists(checkpointFile))
    backupCheckpointFile(checkpointFile);
  boost::filesystem::ofstream ofs(checkpointFile);
  boost::archive::text_oarchive ar(ofs);
  std::cout << "Saving checkpoint to    : " << checkpointFile << '\n';
  boost::serialization::serializeSDPSolverState(ar, x, X, y, Y);
  timers["Last checkpoint"].start();
}

void SDP_Solver::loadCheckpoint(const boost::filesystem::path &checkpointFile)
{
  boost::filesystem::ifstream ifs(checkpointFile);
  boost::archive::text_iarchive ar(ifs);
  std::cout << "Loading checkpoint from : " << checkpointFile << '\n';
  boost::serialization::serializeSDPSolverState(ar, x, X, y, Y);
}

void SDP_Solver::saveSolution(const SDP_Solver_Terminate_Reason terminateReason,
                              const boost::filesystem::path &outFile)
{
  boost::filesystem::ofstream ofs(outFile);
  std::cout << "Saving solution to      : " << outFile << '\n';
  ofs.precision(static_cast<int>(primalObjective.get_prec() * 0.31 + 5));
  ofs << "terminateReason = \"" << terminateReason << "\";\n";
  ofs << "primalObjective = " << primalObjective << ";\n";
  ofs << "dualObjective   = " << dualObjective << ";\n";
  ofs << "dualityGap      = " << dualityGap << ";\n";
  ofs << "primalError     = " << primalError << ";\n";
  ofs << "dualError       = " << dualError << ";\n";
  ofs << "y = " << y << ";\n";
  // ofs << "Y = " << Y << ";\n";
  ofs << "x = " << x << ";\n";
  // ofs << "X = " << X << ";\n";
}
