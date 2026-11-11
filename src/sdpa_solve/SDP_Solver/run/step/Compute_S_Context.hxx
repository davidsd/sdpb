#pragma once

#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpa_solve/Block_Info.hxx"
#include "sdpa_solve/SDP_Solver/run/step/Initialize_P_Context.hxx"
#include "Solver_Run_Config.hxx"

namespace Sdpb::Sdpa
{
  struct Compute_S_Context
  {
    Initialize_P_Context initialize_P_context;
    BigInt_Shared_Memory_Syrk_Context syrk_S_context;
    // TODO
    // Compute_S_Context(const Initialize_P_Config &initialize_p_config,
    //                   const Bigint_Syrk_Config &syrk_S_config)
    //     : initialize_P_context(initialize_p_config),
    //       syrk_S_context(syrk_S_config)
    // {}
    static Compute_S_Context
    create(const Block_Info &block_info, const Solver_Run_Config &cfg,
           const Verbosity verbosity)
    {
      // TODO: reuse shared buffers for trmm and syrk?
      return {Initialize_P_Context(cfg.initialize_P_config),
              BigInt_Shared_Memory_Syrk_Context(
                cfg.syrk_S_config, block_info.node_group_index(),
                block_info.block_indices, verbosity)};
    }
  };
}
