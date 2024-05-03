#pragma once

#include "sdp_solve/Block_Diagonal_Matrix.hxx"
#include "sdpb_util/Timers/Timers.hxx"

// R = mu * I - XY
// R_error = maxAbs(R)
[[nodiscard]] inline El::BigFloat
compute_R_error(const El::BigFloat &mu, const Block_Diagonal_Matrix &minus_XY,
                Timers &timers)
{
  Scoped_Timer R_error_timer(timers, "R_error");
  El::BigFloat R_error = 0;

  for(const auto &block : minus_XY.blocks)
    {
      for(int iLoc = 0; iLoc < block.LocalHeight(); ++iLoc)
        for(int jLoc = 0; jLoc < block.LocalWidth(); ++jLoc)
          {
            const int i = block.GlobalRow(iLoc);
            const int j = block.GlobalCol(jLoc);
            // -XY + mu*I
            auto R_element_value
              = block.GetLocal(iLoc, jLoc) + (i == j ? mu : 0);
            R_error = std::max(R_error, El::Abs(R_element_value));
          }
    }
  return El::mpi::AllReduce(R_error, El::mpi::MAX, El::mpi::COMM_WORLD);
}
