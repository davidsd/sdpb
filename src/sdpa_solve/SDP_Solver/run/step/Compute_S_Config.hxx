#pragma once

#include "sdp_solve/SDP_Solver/run/bigint_syrk/Bigint_Syrk_Config.hxx"
#include "Initialize_P_Config.hxx"
#include "sdp_solve/Solver_Parameters.hxx"
#include "sdpa_solve/Block_Info.hxx"

namespace Sdpb::Sdpa
{
  struct Compute_S_Config
  {
    Initialize_P_Config initialize_P_config;
    Bigint_Syrk_Config syrk_S_config;
    // Total memory estimates

    // Local memory allocated on a node
    [[nodiscard]] size_t node_local_bytes() const;
    // Local memory allocated on a node
    [[nodiscard]] size_t node_shmem_bytes() const;
    // All memory allocated on a node
    [[nodiscard]] size_t node_total_bytes() const;
  };

  Compute_S_Config
  get_compute_S_config(const Environment &env, const Block_Info &block_info,
                       const Solver_Parameters &parameters, size_t other_mem);
}
