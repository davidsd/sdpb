#include "../Outer_Parameters.hxx"

boost::property_tree::ptree to_property_tree(const Outer_Parameters &p)
{
  boost::property_tree::ptree result(to_property_tree(p.solver));

  result.put("sdp", p.functions_path.string());
  result.put("outDir", p.out_directory.string());
  result.put("noFinalCheckpoint", p.no_final_checkpoint);
  result.put("dualityGapReduction", p.duality_gap_reduction);
  result.put("writeSolution", p.write_solution);
  result.put("verbosity", static_cast<int>(p.verbosity));

  return result;
}
