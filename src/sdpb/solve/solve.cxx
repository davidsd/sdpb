//=======================================================================
// Copyright 2014-2015 David Simmons-Duffin.
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
//  http://opensource.org/licenses/MIT)
//=======================================================================

#include "SDP_Solver.hxx"
#include "../Timers.hxx"

#include <El.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>

Timers
solve(const boost::filesystem::path &sdp_directory,
      const boost::filesystem::path &out_file,
      const boost::filesystem::path &checkpoint_in,
      const boost::filesystem::path &checkpoint_out,
      const Block_Info &block_info, const SDP_Solver_Parameters &parameters)
{
  if(El::mpi::Rank() == 0)
    {
      std::cout << "SDPB started at "
                << boost::posix_time::second_clock::local_time() << '\n';
      std::cout << "SDP directory   : " << sdp_directory << '\n';
      std::cout << "out file        : " << out_file << '\n'
                << "checkpoint in   : " << checkpoint_in << '\n'
                << "checkpoint out  : " << checkpoint_out << '\n'
                << "\nParameters:\n"
                << parameters << '\n';
    }
  // Read an SDP from sdpFile and create a solver for it
  El::Grid grid(block_info.mpi_comm);
  SDP sdp(sdp_directory, block_info, grid);
  SDP_Solver solver(parameters, block_info, grid,
                    sdp.dual_objective_b.Height());

  if(exists(checkpoint_in))
    {
      solver.load_checkpoint(checkpoint_in);
    }

  Timers timers(parameters.debug);
  SDP_Solver_Terminate_Reason reason = solver.run(
    parameters, checkpoint_out, block_info, sdp, grid, timers);

  if(El::mpi::Rank() == 0)
    {
      std::cout.precision(
        std::ceil(El::gmp::Precision() * std::log(2.0) / std::log(10.0)));
      std::cout << "-----" << reason << "-----\n"
                << '\n'
                << "primalObjective = " << solver.primal_objective << '\n'
                << "dualObjective   = " << solver.dual_objective << '\n'
                << "dualityGap      = " << solver.duality_gap << '\n'
                << "primalError     = " << solver.primal_error << '\n'
                << "dualError       = " << solver.dual_error << '\n'
                << '\n';
    }

  if(!parameters.no_final_checkpoint)
    {
      solver.save_checkpoint(checkpoint_out);
    }
  solver.save_solution(reason, out_file);
  return timers;
}
