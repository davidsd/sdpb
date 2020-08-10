//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "../sdp_solve.hxx"
#include "../set_stream_precision.hxx"

#include <El.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>

Timers
solve(const Block_Info &block_info, const SDP_Solver_Parameters &parameters)
{
  // Read an SDP from sdpFile and create a solver for it
  El::Grid grid(block_info.mpi_comm.value);
  SDP sdp(parameters.sdp_directory, block_info, grid);
  SDP_Solver solver(parameters, block_info, grid,
                    sdp.dual_objective_b.Height());

  Timers timers(parameters.verbosity >= Verbosity::debug);
  SDP_Solver_Terminate_Reason reason
    = solver.run(parameters, block_info, sdp, grid, timers);

  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      set_stream_precision(std::cout);
      std::cout << "-----" << reason << "-----\n"
                << '\n'
                << "primalObjective = " << solver.primal_objective << '\n'
                << "dualObjective   = " << solver.dual_objective << '\n'
                << "dualityGap      = " << solver.duality_gap << '\n'
                << "primalError     = " << solver.primal_error() << '\n'
                << "dualError       = " << solver.dual_error << '\n'
                << '\n';
    }

  if(!parameters.no_final_checkpoint)
    {
      solver.save_checkpoint(parameters);
    }
  solver.save_solution(reason, timers.front(), parameters.out_directory,
                       parameters.write_solution,
                       block_info.block_indices,
                       parameters.verbosity);
  return timers;
}
