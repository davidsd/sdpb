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

// FIXME: Pass this around instead of having a global.
Timers timers;

int solve(const std::vector<boost::filesystem::path> &sdp_files,
          const boost::filesystem::path &out_file,
          const boost::filesystem::path &checkpoint_file_in,
          const boost::filesystem::path &checkpoint_file_out,
          const SDP_Solver_Parameters &parameters)
{
  if(El::mpi::Rank() == 0)
    {
      std::cout << "SDPB started at "
                << boost::posix_time::second_clock::local_time() << '\n';
      for(auto const &sdpFile : sdp_files)
        {
          std::cout << "SDP file        : " << sdpFile << '\n';
        }
      std::cout << "out file        : " << out_file << '\n'
                << "checkpoint in   : " << checkpoint_file_in << '\n'
                << "checkpoint out  : " << checkpoint_file_out << '\n'
                << "\nParameters:\n"
                << parameters << '\n';
    }
  // Read an SDP from sdpFile and create a solver for it
  SDP_Solver solver(sdp_files, parameters);

  if(exists(checkpoint_file_in))
    {
      solver.load_checkpoint(checkpoint_file_in);
    }

  timers["Solver runtime"].start();
  timers["Last checkpoint"].start();
  SDP_Solver_Terminate_Reason reason = solver.run(checkpoint_file_out);
  timers["Solver runtime"].stop();

  if(El::mpi::Rank() == 0)
    {
      std::cout << "-----" << std::setfill('-') << std::setw(116) << std::left
                << reason << '\n'
                << '\n'
                << "primalObjective = " << solver.primal_objective_elemental
                << '\n'
                << "dualObjective   = " << solver.dual_objective_elemental
                << '\n'
                << "dualityGap      = " << solver.duality_gap_elemental << '\n'
                << "primalError     = " << solver.primal_error_elemental
                << '\n'
                << "dualError       = " << solver.dual_error_elemental << '\n'
                << '\n';
    }

  if(!parameters.no_final_checkpoint)
    {
      solver.save_checkpoint(checkpoint_file_out);
    }
  timers["Last checkpoint"].stop();
  solver.save_solution(reason, out_file);

  if(El::mpi::Rank() == 0)
    {
      std::cout << '\n' << timers;
    }
  timers.write_profile(out_file.string() + ".profiling");

  return 0;
}
