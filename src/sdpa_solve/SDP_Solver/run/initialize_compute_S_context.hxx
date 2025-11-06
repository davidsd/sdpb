#pragma once

#include "sdp_solve/SDP_Solver/run/get_syrk_Q_config.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpa_solve/Block_Info.hxx"
#include "sdpa_solve/SDP_Solver/run/step/Initialize_P_Context.hxx"
#include "step/Compute_S_Config.hxx"
#include "step/Initialize_P_Config.hxx"

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
  };

  inline Compute_S_Context
  initialize_compute_S_context(const Block_Info &block_info,
                               const Compute_S_Config &cfg,
                               const Verbosity verbosity)
  {
    const auto group_index = block_info.node_group_index();
    std::vector<int> group_comm_sizes;
    std::vector<int> blocks_height_per_group;
    for(const auto &group :
        block_info.block_mapping.mapping.at(block_info.node_index))
      {
        group_comm_sizes.push_back(group.num_procs);

        size_t height = 0;
        for(const auto &block_index : group.block_indices)
          {
            const auto dim = block_info.block_dimensions.at(block_index);
            height += dim * dim;
          }
        blocks_height_per_group.push_back(height);
      }

    // TODO: reuse shared buffers for trmm and syrk?
    return {Initialize_P_Context(cfg.initialize_P_config),
            BigInt_Shared_Memory_Syrk_Context(cfg.syrk_S_config, group_index,
                                              block_info.block_indices,
                                              verbosity)};
  }
}
