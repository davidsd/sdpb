#include "Outer_Parameters.hxx"
#include "Function.hxx"
#include "sdp_read/sdp_read.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "sdpb_util/Environment.hxx"

#include "sdpb_util/ostream/ostream_vector.hxx"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace fs = std::filesystem;

// We convert the optimization problem into a regular linear
// programming problem.
//
// 1) Each polynomial in Block::polys adds two variables.  The
// weights of those polynomials are unbounded, while linear
// programming requires a strictly positive variable.  So we
// substitute W_n = w_n+ - w_n+, where both w_n+ and w_n- are
// strictly positive.  All blocks should have the same number of
// polynomials, so we only need to look at the first one.
//
// 2) Each constraint adds one 'slack' variable s_n.  There is one
// constraint per coordinate 'x', with each block having multiple,
// separate coordinates.
//
// 3) One more global variable delta, which gives the linear
// program something to minimize.
//
// This turns the problem
//
//   A_0 . W > 0
//   A_1 . W > 0
//   ...
//
// into
//
//   min delta
//
// where
//
//   A_0 . (w_+ - w_-) + delta - s_0 = 0
//   A_1 . (w_+ - w_-) + delta - s_1 = 0
//   ...
//   w_n+, w_n-, s_n, delta >= 0
//
// There is a constraint for every point that is sampled.  At the
// beginning, we sample the min and max for each block, so there are
// 2*num_blocks constraints.  For the single correlator example, there
// is an additional constraint on the first block at x=0.  In general
//
// => num_rows == num_constraints
// => num_columns == 2*num_weights + num_constraints + 1

std::vector<El::BigFloat> load_vector(const fs::path &vector_path);

void read_function_blocks(
  const fs::path &input_file, std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization,
  std::vector<std::vector<std::vector<std::vector<Function>>>> &functions);

void read_points(const fs::path &input_path,
                 std::vector<std::vector<El::BigFloat>> &points);

std::vector<El::BigFloat> compute_optimal(
  const std::vector<std::vector<std::vector<std::vector<Function>>>> &functions,
  const std::vector<std::vector<El::BigFloat>> &initial_points,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization, const Environment &env,
  const Outer_Parameters &parameters_in,
  const std::chrono::time_point<std::chrono::high_resolution_clock>
    &start_time);

int main(int argc, char **argv)
{
  Environment env(argc, argv);
  Outer_Parameters parameters(argc, argv);
  if(!parameters.is_valid())
    {
      return 0;
    }

  const int64_t precision(parameters.solver.precision);
  El::gmp::SetPrecision(precision);
  // El::gmp wants base-2 bits, but boost::multiprecision wants
  // base-10 digits.
  Boost_Float::default_precision(precision * log(2) / log(10));
  auto start_time = std::chrono::high_resolution_clock::now();
  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      std::cout << boost::posix_time::second_clock::local_time()
                << " Start Outer_Limits" << '\n'
                << "Version: " << SDPB_VERSION_STRING << '\n'
                << parameters << std::endl;
    }

  std::vector<El::BigFloat> objectives, normalization;
  std::vector<std::vector<std::vector<std::vector<Function>>>> functions;
  read_function_blocks(parameters.functions_path, objectives, normalization,
                       functions);

  std::vector<std::vector<El::BigFloat>> initial_points;
  read_points(parameters.points_path, initial_points);

  std::vector<El::BigFloat> weights(
    compute_optimal(functions, initial_points, objectives, normalization, env,
                    parameters, start_time));

  El::BigFloat optimal(0);
  for(size_t index(0); index < objectives.size(); ++index)
    {
      optimal += objectives[index] * weights[index];
    }
  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      set_stream_precision(std::cout);
      std::cout << "optimal: " << optimal << "\n";
    }
  if(El::mpi::Rank() == 0)
    {
      if(parameters.verbosity >= Verbosity::regular)
        {
          std::cout << "Saving solution to " << parameters.output_path << "\n";
        }
      fs::create_directories(parameters.output_path.parent_path());
      std::ofstream output(parameters.output_path);
      set_stream_precision(output);
      output << "{\n  \"optimal\": \"" << optimal << "\",\n"
             << "  \"y\":\n  [\n";
      for(auto weight(weights.begin()); weight != weights.end(); ++weight)
        {
          if(weight != weights.begin())
            {
              output << ",\n";
            }
          output << "    \"" << *weight << "\"";
        }

      output << "\n  ],\n  \"options\": \n";
      boost::property_tree::write_json(output, to_property_tree(parameters));
      output << "}\n";
    }
}
