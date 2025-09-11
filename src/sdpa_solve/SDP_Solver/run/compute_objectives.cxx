#include "sdpa_solve/SDP.hxx"
#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/Timers/Timers.hxx"

namespace Sdpb::Sdpa
{
  void compute_objectives(const SDP &sdp, const Primal_Dist_Vector &x,
                          const Block_Diagonal_Matrix &Y,
                          El::BigFloat &primal_objective,
                          El::BigFloat &dual_objective,
                          El::BigFloat &duality_gap, Timers &timers)
  {
    Scoped_Timer objectives_timer(timers, "objectives");
    // P-obj = c.x
    primal_objective = El::Dotu(sdp.primal_objective_c, x);
    // D-obj = Tr(F_0 * Y)
    dual_objective = dotu(sdp.sdp_block_F_0, Y);
    El::mpi::Broadcast(dual_objective, 0, El::mpi::COMM_WORLD);

    duality_gap
      = Abs(primal_objective - dual_objective)
        / Max(Abs(primal_objective) + Abs(dual_objective), El::BigFloat(1));
  }
}
