#include "Function.hxx"
#include "../sdp_read.hxx"
#include "../sdp_solve.hxx"

#include "../ostream_vector.hxx"

#include <boost/date_time/posix_time/posix_time.hpp>

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

std::vector<El::BigFloat>
load_vector(const boost::filesystem::path &vector_path);

void read_function_blocks(
  const boost::filesystem::path &input_file,
  std::vector<El::BigFloat> &objectives,
  std::vector<El::BigFloat> &normalization,
  std::vector<std::vector<std::vector<std::vector<Function>>>> &functions);

void read_points(const boost::filesystem::path &input_path,
                 std::vector<std::vector<El::BigFloat>> &points);

void convert_matrices_to_functions(
  const El::BigFloat &max_delta,
  const std::vector<Positive_Matrix_With_Prefactor> &matrices,
  std::vector<std::vector<std::vector<std::vector<Function>>>> &functions);

std::vector<El::BigFloat> compute_optimal(
  const std::vector<std::vector<std::vector<std::vector<Function>>>> &functions,
  const std::vector<std::vector<El::BigFloat>> &initial_points,
  const std::vector<El::BigFloat> &objectives,
  const std::vector<El::BigFloat> &normalization,
  const Solver_Parameters &parameters_in);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);
  Solver_Parameters parameters(argc, argv);
  if(!parameters.is_valid())
    {
      return 0;
    }

  const int64_t precision(parameters.precision);
  El::gmp::SetPrecision(precision);
  // El::gmp wants base-2 bits, but boost::multiprecision wants
  // base-10 digits.
  Boost_Float::default_precision(precision * log(2) / log(10));

  if(parameters.verbosity >= Verbosity::regular && El::mpi::Rank() == 0)
    {
      std::cout << "Outer_Limits started at "
                << boost::posix_time::second_clock::local_time() << '\n'
                << parameters << '\n'
                << std::flush;
    }

  {
    std::vector<El::BigFloat> objectives, normalization;
    std::vector<Positive_Matrix_With_Prefactor> matrices;
    std::vector<std::vector<El::BigFloat>> initial_points;

    // read_input("outer_small.nsv", objectives, normalization, matrices);
    // read_points("outer_small_points.json", initial_points);
    // // read_points("outer_small_points_extra.json", initial_points);

    read_input("outer_2x2.nsv", objectives, normalization, matrices);
    read_points("outer_2x2_points.json", initial_points);

    // read_input("test/spectrum_test.json", objectives, normalization,
    // matrices); read_input("test/toy_damped_duplicate.json", objectives,
    // normalization, matrices);
    // read_input("test/toy_damped_3.json", objectives, normalization,
    // matrices);

    // read_input("test/toy_damped.json", objectives, normalization, matrices);
    // read_points("test/toy_damped_points.json", initial_points);

    std::vector<std::vector<std::vector<std::vector<Function>>>> functions;
    const El::BigFloat max_delta(64); // This is completely arbitrary.
    convert_matrices_to_functions(max_delta, matrices, functions);

    std::vector<El::BigFloat> weights(compute_optimal(
      functions, initial_points, objectives, normalization, parameters));

    El::BigFloat optimal(0);
    for(size_t index(0); index < objectives.size(); ++index)
      {
        optimal += objectives[index] * weights[index];
      }
    if(El::mpi::Rank() == 0)
      {
        std::cout.precision(precision / 3.3);
        std::cout << "optimal: " << optimal << " " << weights << "\n";
      }
  }
}
