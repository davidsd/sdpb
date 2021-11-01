#include "../Dynamic_Parameters.hxx"

boost::property_tree::ptree to_property_tree(const Dynamic_Parameters &p)
{
  boost::property_tree::ptree result(to_property_tree(p.solver));

  result.put("sdpDir", p.sdp_path.string());
  result.put("newsdpDir", p.new_sdp_path.string());
  result.put("outDir", p.out_directory.string());
  result.put("noFinalCheckpoint", p.no_final_checkpoint);
  result.put("writeSolution", p.write_solution);
  result.put("procsPerNode", p.procs_per_node);
  result.put("procGranularity", p.proc_granularity);
  result.put("verbosity", static_cast<int>(p.verbosity));
  result.put("StepSizeAlpha", p.alpha);
  result.put("NumberExternalParameters", p.N_external_parameters);

  return result;
}
