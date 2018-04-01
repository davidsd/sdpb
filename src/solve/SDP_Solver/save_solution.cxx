#include "../SDP_Solver.hxx"

void SDP_Solver::save_solution(const SDP_Solver_Terminate_Reason terminate_reason,
                              const boost::filesystem::path &out_file)
{
  boost::filesystem::ofstream ofs(out_file);
  std::cout << "Saving solution to      : " << out_file << '\n';
  ofs.precision(static_cast<int>(primal_objective.get_prec() * 0.31 + 5));
  ofs << "terminateReason = \"" << terminate_reason << "\";\n";
  ofs << "primalObjective = " << primal_objective << ";\n";
  ofs << "dualObjective   = " << dual_objective << ";\n";
  ofs << "dualityGap      = " << duality_gap << ";\n";
  ofs << "primalError     = " << primal_error << ";\n";
  ofs << "dualError       = " << dual_error << ";\n";
  ofs << "y = " << y << ";\n";
  // ofs << "Y = " << Y << ";\n";
  ofs << "x = " << x << ";\n";
  // ofs << "X = " << X << ";\n";
}
