#include "Output_SDP/Output_SDP.hxx"
#include "Pmp2sdp_Parameters/Pmp2sdp_Parameters.hxx"
#include "Dual_Constraint_Group.hxx"
#include "write_sdp.hxx"
#include "pmp_read/pmp_read.hxx"
#include "apply_encoding_points.hxx"

void set_block_data_output_precision(mp_bitcnt_t precision);
#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/Timers/Timers.hxx"
#include "sdpb_util/assert.hxx"

#include <string>
#include <boost/program_options.hpp>
#include <filesystem>

namespace fs = std::filesystem;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      Pmp2sdp_Parameters parameters(argc, argv);
      if(!parameters.is_valid())
        {
          return 0;
        }
      // Print command line
      if(parameters.verbosity >= Verbosity::debug && El::mpi::Rank() == 0)
        {
          std::vector<std::string> arg_list(argv, argv + argc);
          for(const auto &arg : arg_list)
            std::cout << arg << " ";
          std::cout << std::endl;
        }

      // --encodingPoints: the whole conversion runs at a higher precision
      // (the step to the encoding points cancels digits), and the SDP is
      // truncated to --precision before it is written.
      const mp_bitcnt_t working_precision
        = parameters.precision
          + (parameters.two_set_encoding()
               ? parameters.encoding_extra_precision
               : 0);
      Environment::set_precision(working_precision);

      Timers timers(env, parameters.verbosity);
      Scoped_Timer timer(timers, "pmp2sdp");

      auto pmp = read_polynomial_matrix_program(env, parameters.input_file,
                                                parameters.max_num_poles,
                                                parameters.verbosity, timers);

      Encoding_Points_Report encoding_report;
      if(parameters.two_set_encoding())
        {
          Scoped_Timer encoding_timer(timers, "encoding_points");
          encoding_report = apply_encoding_points(pmp, parameters);
        }

      Output_SDP sdp(pmp, parameters.command_arguments, timers);
      if(parameters.two_set_encoding())
        {
          double digits_lost = 0;
          for(const auto &group : sdp.dual_constraint_groups)
            digits_lost = std::max(digits_lost, group.encoding_digits_lost);
          digits_lost = El::mpi::AllReduce(digits_lost, El::mpi::MAX,
                                           El::mpi::COMM_WORLD);
          const double reach = El::mpi::AllReduce(
            encoding_report.max_reach, El::mpi::MAX, El::mpi::COMM_WORLD);
          const double extra_digits
            = parameters.encoding_extra_precision * std::log10(2.0);
          if(El::mpi::Rank() == 0
             && parameters.verbosity >= Verbosity::regular)
            El::Output("encodingPoints=", parameters.encoding_points,
                       ": largest encoding point ", reach,
                       " times the largest approximation point; at most ",
                       digits_lost,
                       " decimal digits lost to cancellation in the step to "
                       "the encoding points (encodingExtraPrecision has ",
                       extra_digits, ")");
          // Leave at least 32 bits (9.6 digits) of the extra precision.
          ASSERT(digits_lost <= extra_digits - 32 * std::log10(2.0),
                 "--encodingPoints: the step to the encoding points lost ",
                 digits_lost, " decimal digits, more than "
                 "--encodingExtraPrecision=",
                 parameters.encoding_extra_precision,
                 " bits minus 32 allow. Raise --encodingExtraPrecision.");

          // Blocks of degree 0 are not re-encoded: bring them to the
          // re-encoded blocks' scale (an equivalent constraint; see
          // Dual_Constraint_Group/scale_degree_0_blocks.cxx).
          for(size_t i = 0; i < pmp.matrices.size(); ++i)
            ASSERT_EQUAL(sdp.dual_constraint_groups.at(i).block_index,
                         pmp.matrix_index_local_to_global.at(i));
          const auto degree_0
            = scale_degree_0_blocks(sdp.dual_constraint_groups, pmp.matrices);
          if(El::mpi::Rank() == 0
             && parameters.verbosity >= Verbosity::regular
             && degree_0.num_degree_0 > 0)
            {
              if(degree_0.num_reencoded == 0)
                El::Output("encodingPoints: ", degree_0.num_degree_0,
                           " blocks of degree 0 and no re-encoded block: "
                           "the degree-0 blocks keep their scale");
              else
                El::Output(
                  "encodingPoints: ", degree_0.num_degree_0,
                  " blocks of degree 0 brought to the median scale of the ",
                  degree_0.num_reencoded,
                  " re-encoded blocks (largest |entry| of B and c 10^",
                  degree_0.median_log2 * std::log10(2.0), "): ",
                  degree_0.num_scaled, " multiplied by 2^e, e from ",
                  degree_0.min_exponent, " to ", degree_0.max_exponent);
            }

          // Write at --precision (binary blocks must have SDPB's precision).
          // El's global precision can be set only once, so each written
          // number is truncated in place and the block files carry
          // --precision.
          const mp_bitcnt_t precision = parameters.precision;
          set_block_data_output_precision(precision);
          truncate_precision(sdp.objective_const, precision);
          truncate_precision(sdp.dual_objective_b, precision);
          if(sdp.normalization.has_value())
            truncate_precision(sdp.normalization.value(), precision);
          for(auto &group : sdp.dual_constraint_groups)
            {
              truncate_precision(group.constraint_matrix, precision);
              truncate_precision(group.constraint_constants, precision);
              for(auto &basis : group.bilinear_bases)
                truncate_precision(basis, precision);
            }
          for(auto &matrix : pmp.matrices)
            {
              truncate_precision(matrix.sample_points, precision);
              truncate_precision(matrix.sample_scalings, precision);
              truncate_precision(matrix.reduced_sample_scalings, precision);
              truncate_precision(matrix.encoding_points, precision);
              truncate_precision(matrix.encoding_scalings, precision);
            }
        }
      write_sdp(parameters.output_path, sdp, pmp, parameters.output_format,
                parameters.zip, timers, parameters.verbosity);
      if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
        {
          El::Output("Processed ", sdp.num_blocks, " SDP blocks in ",
                     (double)timer.timer().elapsed_milliseconds() / 1000,
                     " seconds, output: ", parameters.output_path.string());
        }

      if(parameters.verbosity >= Verbosity::debug)
        {
          timers.write_profile(parameters.output_path.string()
                               + ".profiling/profiling."
                               + std::to_string(El::mpi::Rank()));
        }
    }
  catch(std::exception &e)
    {
      std::cerr << "Error: " << e.what() << "\n" << std::flush;
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
  catch(...)
    {
      std::cerr << "Unknown Error\n" << std::flush;
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
