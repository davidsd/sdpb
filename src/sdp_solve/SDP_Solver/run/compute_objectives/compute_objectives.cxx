#include "../../../SDP.hxx"
#include "../../../../Timers.hxx"

El::BigFloat dot(const Block_Vector &a, const Block_Vector &b);

void compute_objectives(const SDP &sdp, const Block_Vector &x,
                        const Block_Vector &y, El::BigFloat &primal_objective,
                        El::BigFloat &dual_objective,
                        El::BigFloat &duality_gap, Timers &timers)
{
  Scoped_Timer objectives_timer(timers, "run.objectives");
  primal_objective = sdp.objective_const + dot(sdp.primal_objective_c, x);
  // dual_objective_b is duplicated amongst the processors.  y is
  // duplicated amongst the blocks, but it is possible for some
  // processors to have no blocks.  In principle, we only need to
  // compute the dot product on the first block, but then we would
  // have to make sure that we compute that product over all
  // processors that own that block.
  if(!y.blocks.empty())
    {
      dual_objective = sdp.objective_const
                       + El::Dotu(sdp.dual_objective_b, y.blocks.front());
    }
  El::mpi::Broadcast(dual_objective, 0, El::mpi::COMM_WORLD);

  duality_gap
    = Abs(primal_objective - dual_objective)
      / Max(Abs(primal_objective) + Abs(dual_objective), El::BigFloat(1));
}
