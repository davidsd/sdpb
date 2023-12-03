#include "approx_objective/Approx_Objective.hxx"

void compute_dx_dy(const Block_Info &block_info, const SDP &d_sdp,
                   const Block_Vector &x, const Block_Vector &y,
                   const Block_Diagonal_Matrix &schur_complement_cholesky,
                   const Block_Matrix &schur_off_diagonal,
                   const El::DistMatrix<El::BigFloat> &Q, Block_Vector &dx,
                   Block_Vector &dy);

// Only linear term
Approx_Objective::Approx_Objective(const SDP &sdp, const SDP &d_sdp,
                                   const Block_Vector &x,
                                   const Block_Vector &y)
    : dd_objective(0)
{
  // b.y
  objective
    = El::Dot(sdp.dual_objective_b, y.blocks.at(0)) + sdp.objective_const;

  El::BigFloat linear(0);
  // dconst
  d_objective = d_sdp.objective_const;
  // db.y
  d_objective += El::Dot(d_sdp.dual_objective_b, y.blocks.at(0));

  El::BigFloat local_linear(0);
  for(size_t block(0); block != x.blocks.size(); ++block)
    {
      // dc.x
      local_linear
        += Dotu(d_sdp.primal_objective_c.blocks.at(block), x.blocks.at(block));

      {
        // temp = dB.y
        El::DistMatrix<El::BigFloat> temp(x.blocks.at(block));
        El::Zero(temp);
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                 d_sdp.free_var_matrix.blocks[block], y.blocks.at(0),
                 El::BigFloat(0.0), temp);

        // -x.dB.y/
        local_linear -= El::Dotu(temp, x.blocks.at(block));
      }
    }
  if(!x.blocks.empty() && x.blocks.at(0).Grid().Rank() != 0)
    {
      local_linear = 0;
    }

  d_objective
    += El::mpi::AllReduce(local_linear, El::mpi::SUM, El::mpi::COMM_WORLD);
  objective += d_objective + dd_objective;
}

// Linear and quadratic terms
Approx_Objective::Approx_Objective(
  const Block_Info &block_info, const SDP &sdp, const SDP &d_sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q)
{
  Block_Vector dx(x), dy(y);
  compute_dx_dy(block_info, d_sdp, x, y, schur_complement_cholesky,
                schur_off_diagonal, Q, dx, dy);

  // b.y
  objective
    = El::Dot(sdp.dual_objective_b, y.blocks.at(0)) + sdp.objective_const;

  El::BigFloat linear(0), quad(0);
  // dconst
  d_objective = d_sdp.objective_const;
  // db.y
  d_objective += El::Dot(d_sdp.dual_objective_b, y.blocks.at(0));

  // db.dy/2
  dd_objective = El::Dot(d_sdp.dual_objective_b, dy.blocks.at(0)) / 2;

  El::BigFloat local_linear(0), local_quad(0);
  for(size_t block(0); block != x.blocks.size(); ++block)
    {
      // dc.x
      local_linear
        += Dotu(d_sdp.primal_objective_c.blocks.at(block), x.blocks.at(block));

      // dc.dx/2
      local_quad
        += Dotu(d_sdp.primal_objective_c.blocks.at(block), dx.blocks.at(block))
           / 2;
      {
        // temp = dB.y
        El::DistMatrix<El::BigFloat> temp(x.blocks.at(block));
        El::Zero(temp);
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                 d_sdp.free_var_matrix.blocks[block], y.blocks.at(0),
                 El::BigFloat(0.0), temp);

        // -x.dB.y/
        local_linear -= El::Dotu(temp, x.blocks.at(block));

        // -dx.dB.y/2
        local_quad -= El::Dotu(temp, dx.blocks.at(block)) / 2;

        // temp = dB.dy
        El::Gemv(El::Orientation::NORMAL, El::BigFloat(1.0),
                 d_sdp.free_var_matrix.blocks[block], dy.blocks.at(0),
                 El::BigFloat(0.0), temp);

        // -x.dB.dy/2
        local_quad -= El::Dotu(temp, x.blocks.at(block)) / 2;
      }
    }
  if(!x.blocks.empty() && x.blocks.at(0).Grid().Rank() != 0)
    {
      local_linear = 0;
      local_quad = 0;
    }

  d_objective
    += El::mpi::AllReduce(local_linear, El::mpi::SUM, El::mpi::COMM_WORLD);
  dd_objective
    += El::mpi::AllReduce(local_quad, El::mpi::SUM, El::mpi::COMM_WORLD);
  objective += d_objective + dd_objective;
}
