#include "pmp_read/pmp_read.hxx"
#include "sdp_solve/read_text_block.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <boost/program_options.hpp>

#include <filesystem>
#include <ostream>
#include <string>

namespace fs = std::filesystem;
namespace po = boost::program_options;

template <class T> struct Vector : std::vector<T>
{
  // Required for El::Matrix<Vector<T>>
  void Zero() { this->clear(); }
};

struct Vector_Diff
{
  size_t index = 0;
  El::BigFloat x = 0;
  El::BigFloat max_diff = 0;
  void Zero()
  {
    index = 0;
    max_diff = 0;
  }
  friend std::ostream &operator<<(std::ostream &os, const Vector_Diff &obj)
  {
    if(obj.max_diff == El::BigFloat(0))
      return os << "0";

    return os << "(max=" << obj.max_diff << ", x=" << obj.x << ")";
  }
};

struct Matrix_Diff
{
  El::Matrix<Vector_Diff> diffs;
  friend std::ostream &operator<<(std::ostream &os, const Matrix_Diff &obj)
  {
    const auto &matrix = obj.diffs;
    const auto label = "";
    os << El::DimsString(matrix, label);
    os << "\n[\n";
    for(int i = 0; i < matrix.Height(); ++i)
      {
        for(int j = 0; j < matrix.Width(); ++j)
          {
            os << " " << matrix.Get(i, j);
          }
        os << "\n";
      }
    os << "]";

    El::BigFloat max_diff = 0;
    for(int i = 0; i < obj.diffs.Height(); ++i)
      for(int j = 0; j < obj.diffs.Width(); ++j)
        {
          max_diff = std::max(max_diff, obj.diffs.CRef(i, j).max_diff);
        }

    os << "\n  max_diff=" << max_diff;
    return os;
  }
};

// For each poly_vec element of PMP matrix,
// compute product poly_vec.z and evaluate at sample points (including prefactor)
[[nodiscard]] El::Matrix<Vector<El::BigFloat>>
sample_pmp_matrix(const Polynomial_Vector_Matrix &matrix,
                  const El::Matrix<El::BigFloat> &z)
{
  El::Matrix<Vector<El::BigFloat>> result(matrix.polynomials.Height(),
                                          matrix.polynomials.Width());
  for(int i = 0; i < matrix.polynomials.Height(); ++i)
    for(int j = 0; j < matrix.polynomials.Width(); ++j)
      {
        const auto &poly_vec = matrix.polynomials.CRef(i, j);

        auto &values = result.Ref(i, j);
        for(int point_index = 0; point_index < matrix.sample_points.size();
            ++point_index)
          {
            const auto &x = matrix.sample_points.at(point_index);
            const auto &scale = matrix.sample_scalings.at(point_index);

            ASSERT_EQUAL(poly_vec.size(), z.Height());
            El::BigFloat value = 0;
            for(int z_index = 0; z_index < z.Height(); ++z_index)
              {
                const auto &poly = poly_vec.at(z_index);
                value += z(z_index, 0) * scale * poly(x);
              }
            values.push_back(value);
          }
      }
  return result;
}

// Attempting to measure difference of two matrices
Matrix_Diff matrix_diff(const El::Matrix<Vector<El::BigFloat>> &left,
                        const El::Matrix<Vector<El::BigFloat>> &right,
                        const std::vector<El::BigFloat> &sample_points)
{
  ASSERT_EQUAL(left.Height(), right.Height());
  ASSERT_EQUAL(left.Width(), right.Width());

  Matrix_Diff result;
  result.diffs.Resize(left.Height(), left.Width());
  for(int i = 0; i < left.Height(); ++i)
    for(int j = 0; j < left.Width(); ++j)
      {
        const auto &left_vec = left.CRef(i, j);
        const auto &right_vec = right.CRef(i, j);
        ASSERT_EQUAL(left_vec.size(), right_vec.size());

        Vector_Diff vec_diff;
        vec_diff.x = sample_points.at(0);

        for(int n = 0; n < left_vec.size(); ++n)
          {
            auto curr_diff = El::Abs(left_vec.at(n) - right_vec.at(n));
            if(curr_diff > vec_diff.max_diff)
              {
                vec_diff.index = n;
                vec_diff.x = sample_points.at(n);
                vec_diff.max_diff = curr_diff;
              }
          }

        result.diffs.Set(i, j, vec_diff);
      }
  return result;
}

Matrix_Diff diff_pvm(const Polynomial_Vector_Matrix &left,
                     const Polynomial_Vector_Matrix &right,
                     const El::Matrix<El::BigFloat> &z)
{
  return matrix_diff(sample_pmp_matrix(left, z), sample_pmp_matrix(right, z),
                     left.sample_points);
}
void compare_pmp(const Polynomial_Matrix_Program &left_pmp,
                 const Polynomial_Matrix_Program &right_pmp,
                 const El::Matrix<El::BigFloat> &z)
{
  const size_t num_matrices = left_pmp.num_matrices;
  std::vector<size_t> left_block_index_global_to_local(num_matrices);
  std::vector<size_t> right_block_index_global_to_local(num_matrices);

  for(size_t index_local = 0; index_local < num_matrices; ++index_local)
    {
      left_block_index_global_to_local.at(
        left_pmp.matrix_index_local_to_global.at(index_local))
        = index_local;
      right_block_index_global_to_local.at(
        right_pmp.matrix_index_local_to_global.at(index_local))
        = index_local;
    }

  for(size_t block_index = 0; block_index < num_matrices; ++block_index)
    {
      const auto left_local_index
        = left_block_index_global_to_local.at(block_index);
      const auto right_local_index
        = right_block_index_global_to_local.at(block_index);

      const auto &left_path = left_pmp.block_paths.at(left_local_index);
      const auto &right_path = right_pmp.block_paths.at(right_local_index);

      const auto &left_pvm = left_pmp.matrices.at(left_local_index);
      const auto &right_pvm = right_pmp.matrices.at(right_local_index);
      auto diff = diff_pvm(left_pvm, right_pvm, z);

      El::Output("\nblock_index=", block_index);
      El::Output("  left=", left_path);
      El::Output("  right=", right_path);
      El::Output("  diff=", diff);
    }
}

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  // There can be problems with different PMP block mapping etc.
  ASSERT_EQUAL(El::mpi::Size(), 1, "TODO: cannot run pmp_compare in parallel");

  try
    {
      int precision;
      fs::path left_pmp_path, right_pmp_path;
      fs::path z_path;
      Verbosity verbosity = Verbosity::trace;

      // Parse command-line arguments
      {
        po::options_description options("Basic options");
        options.add_options()("help,h", "Show this helpful message.");
        options.add_options()("left,l",
                              po::value<fs::path>(&left_pmp_path)->required(),
                              "first PMP");
        options.add_options()("right,r",
                              po::value<fs::path>(&right_pmp_path)->required(),
                              "second PMP");
        options.add_options()(
          "z", po::value<fs::path>(&z_path)->required(),
          "Path to z.txt containing a vector that should be "
          "dotted with each PMP element.");
        options.add_options()(
          "precision", po::value<int>(&precision)->required(),
          "The precision, in the number of bits, for numbers in the "
          "computation. ");
        options.add_options()(
          "verbosity", po::value<Verbosity>(&verbosity),
          "Verbosity (0,1,2,3 or none,regular,debug,trace)");

        po::positional_options_description positional;
        positional.add("precision", 1);
        positional.add("left", 1);
        positional.add("right", 1);
        positional.add("z", 1);

        po::variables_map variables_map;
        po::store(po::command_line_parser(argc, argv)
                    .options(options)
                    .positional(positional)
                    .run(),
                  variables_map);

        if(El::mpi::Rank() == 0 && variables_map.count("help") != 0)
          {
            El::Output("pmp_compare reads two PMPs and vector z.");
            El::Output("Then it multiplies each polynomial vector in PMPs by "
                       "z, and evaluates them at sample points.");
            El::Output("After that, each PMP matrix is a matrix of vectors of "
                       "numbers.");
            El::Output("The program prints difference between each matrix.");
            El::Output(options);
            return 0;
          }

        po::notify(variables_map);

        // Print command line
        if(verbosity >= Verbosity::debug && El::mpi::Rank() == 0)
          {
            std::vector<std::string> arg_list(argv, argv + argc);
            for(const auto &arg : arg_list)
              std::cout << arg << " ";
            std::cout << std::endl;
          }

        ASSERT(fs::exists(left_pmp_path),
               "--left PMP does not exist: ", left_pmp_path);
        ASSERT(fs::exists(right_pmp_path),
               "--right PMP does not exist: ", right_pmp_path);
        ASSERT(fs::exists(z_path), "z.txt does not exist: ", z_path);
      }

      Environment::set_precision(precision);

      Timers timers(env, verbosity);
      Scoped_Timer timer(timers, "pmp_compare");

      auto left_pmp
        = read_polynomial_matrix_program(env, left_pmp_path, timers);
      auto right_pmp
        = read_polynomial_matrix_program(env, right_pmp_path, timers);

      const size_t num_blocks = left_pmp.num_matrices;
      ASSERT_EQUAL(num_blocks, right_pmp.num_matrices);

      El::Matrix<El::BigFloat> z(left_pmp.objective.size(), 1);
      read_text_block(z, z_path);

      Scoped_Timer compare_timer(timers, "compare");
      compare_pmp(left_pmp, right_pmp, z);
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
