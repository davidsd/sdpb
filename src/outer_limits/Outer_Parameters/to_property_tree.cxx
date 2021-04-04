#include "../Outer_Parameters.hxx"

boost::property_tree::ptree to_property_tree(const Outer_Parameters &p)
{
  boost::property_tree::ptree result(to_property_tree(p.solver));

  result.put("functions", p.functions_path.string());
  result.put("points", p.points_path.string());
  result.put("out", p.output_path.string());
  result.put("noFinalCheckpoint", p.no_final_checkpoint);
  result.put("dualityGapReduction", p.duality_gap_reduction);
  result.put("writeSolution", p.write_solution);
  result.put("verbosity", static_cast<int>(p.verbosity));

  return result;
}
