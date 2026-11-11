#pragma once

#include "SDPB_Parameters.hxx"
#include "sdpb_util/Environment.hxx"

#include <filesystem>

inline std::filesystem::path
get_block_timings_path(const SDPB_Parameters &parameters)
{
  std::filesystem::path block_timings_path;
  const auto sdp_block_timings_path = parameters.sdp_path / "block_timings";
  const auto checkpoint_block_timings_path
    = parameters.solver.checkpoint_in / "block_timings";

  if(exists(checkpoint_block_timings_path))
    block_timings_path = checkpoint_block_timings_path;

  if(exists(sdp_block_timings_path))
    block_timings_path = sdp_block_timings_path;
  return block_timings_path;
}

template <class TSolver>
typename TSolver::Block_Info_Type
create_block_info(const Environment &env, const SDPB_Parameters &parameters,
                  const std::filesystem::path &block_timings_path)
{
  return TSolver::Block_Info_Type::create(
    env, parameters.sdp_path, block_timings_path, parameters.proc_granularity,
    parameters.verbosity);
}

template <class TSolver>
typename TSolver::Block_Info_Type
create_block_info(const Environment &env, const SDPB_Parameters &parameters,
                  const El::Matrix<int32_t> &block_timings_ms)
{
  return TSolver::Block_Info_Type::create(
  env, parameters.sdp_path, block_timings_ms, parameters.proc_granularity,
  parameters.verbosity);
}

