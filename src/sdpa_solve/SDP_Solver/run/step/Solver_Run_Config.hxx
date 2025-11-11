#pragma once

#include "sdp_solve/SDP_Solver/run/bigint_syrk/Bigint_Syrk_Config.hxx"
#include "Initialize_P_Config.hxx"
#include "sdp_solve/Solver_Parameters.hxx"
#include "sdpa_solve/Block_Info.hxx"

namespace Sdpb::Sdpa
{
  struct SDP_Solver;
  struct SDP;
  struct Solver_Run_Config
  {
    // Precompute various matrix sizes
    struct Sizes
    {
      size_t X_size;
      size_t X_bytes;
      size_t x_bytes;
      size_t S_size;
      size_t S_bytes;
      size_t SDP_bytes;
      size_t SDP_solver_bytes;
      // Total height of P matrix blocks for each MPI group
      std::vector<size_t> P_group_heights;
      Sizes(const Block_Info &block_info, const SDP &sdp,
            const SDP_Solver &solver);
    };

    Sizes sizes;
    Initialize_P_Config initialize_P_config;
    Bigint_Syrk_Config syrk_S_config;

    Solver_Run_Config(const Sizes &sizes,
                      const Initialize_P_Config &initialize_P_config,
                      const Bigint_Syrk_Config &syrk_S_config);
    static Solver_Run_Config
    create(const Environment &env, const Block_Info &block_info,
           const SDP &sdp, const SDP_Solver &solver,
           const Solver_Parameters &parameters, const Verbosity &verbosity);
  };
}
