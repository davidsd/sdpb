#include "../SDP_Solver.hxx"

void SDP_Solver::save_solution(
  const SDP_Solver_Terminate_Reason terminate_reason,
  const boost::filesystem::path &out_file)
{
  boost::filesystem::ofstream ofs(out_file);
  std::cout << "Saving solution to      : " << out_file << '\n';
  ofs << "terminateReason = \"" << terminate_reason << "\";\n"
      << "primalObjective = " << primal_objective_elemental << ";\n"
      << "dualObjective   = " << dual_objective_elemental << ";\n"
      << "dualityGap      = " << duality_gap_elemental << ";\n"
      << "primalError     = " << primal_error_elemental << ";\n"
      << "dualError       = " << dual_error_elemental << ";\n"
      << "y = {";
  El::Print(y_elemental, "", ofs);
  ofs << "};\nx = {";
  for(auto &block : x_elemental.blocks)
    {
      El::Print(block, "", ofs);
    }
  ofs << "};\n";
}
