#pragma once

#include "sdp_solve/Block_Matrix/Block_Matrix.hxx"
#include "sdpa_solve/Block_Info.hxx"
#include "sdpa_solve/SDP.hxx"
#include "sdpa_solve/SDP_Solver/run/step/Initialize_P_Context.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/initialize_bigint_syrk_context.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/bigint_shared_memory/Vector_Matrix_Residues_Window.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/Bigint_Syrk_Config.hxx"
#include "step/Compute_S_Config.hxx"
#include "step/Initialize_P_Config.hxx"

namespace Sdpb::Sdpa
{
  struct Compute_S_Context
  {
    Initialize_P_Context initialize_P_context;
    BigInt_Shared_Memory_Syrk_Context syrk_P_context;
    // TODO
    // Compute_S_Context(const Initialize_P_Config &initialize_p_config,
    //                   const Bigint_Syrk_Config &syrk_P_config)
    //     : initialize_P_context(initialize_p_config),
    //       syrk_P_context(syrk_P_config)
    // {}
  };

  inline Compute_S_Context
  initialize_compute_S_context(const Block_Info &block_info,
                               const Compute_S_Config &cfg,
                               const Verbosity verbosity)
  {
    const auto group_index = block_info.node_group_index();
    const auto precision = El::gmp::Precision();

    const auto block_width = block_info.primal_dimension;

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
            BigInt_Shared_Memory_Syrk_Context(cfg.syrk_P_config, group_index,
                                              block_info.block_indices,
                                              verbosity)};
  }
}
