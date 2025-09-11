#include "sdpa_solve/SDP_Solver.hxx"

namespace Sdpb::Sdpa
{
  // dual_residues[p] = c[p] - F[p,a,b] Y[a,b]
  void
  compute_dual_residues_and_error(const SDP &sdp,
                                  const Block_Diagonal_Matrix &Y,
                                  Primal_Dist_Vector &dual_residues,
                                  El::BigFloat &dual_error, Timers &timers)
  {
    Scoped_Timer dual_residues_timer(timers, "computeDualResidues");
    dual_residues = sdp.primal_objective_c;
    for(int i = 0; i < dual_residues.Height(); ++i)
      {
        const auto d = dual_residues.Get(i, 0);
        dual_residues.Set(i, 0, d - dotu(sdp.sdp_blocks_F.at(i), Y));
      }
    dual_error = El::MaxAbs(dual_residues);
  }
}
