#include "pmp_read/pmp_read.hxx"
#include "sdp_solve/read_text_block.hxx"
#include "sdpb_util/copy_matrix.hxx"
#include "sdpb_util/Timers/Timers.hxx"
#include "sdpb_util/ostream/ostream_vector.hxx"

#include <filesystem>
#include <ostream>
#include <string>
#include <boost/program_options.hpp>
#include <rapidjson/filewritestream.h>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/rapidjson.h>
#include <rapidjson/writer.h>

namespace fs = std::filesystem;
namespace po = boost::program_options;

template <class T> struct Vector : std::vector<T>
{
  // Required for El::Matrix<Vector<T>>
  void Zero() { this->clear(); }
};

struct Matrix_Diff
{
  El::Matrix<El::BigFloat> max_diff;
  El::Matrix<El::BigFloat> max_x;

  void resize(const int height, const int width)
  {
    max_diff.Resize(height, width);
    El::Zero(max_diff);
    max_x.Resize(height, width);
    El::Zero(max_x);
  }

  [[nodiscard]] std::tuple<int, int, El::BigFloat, El::BigFloat>
  get_max_element_and_x() const
  {
    auto [i, j, value] = El::MaxAbsLoc(max_diff);
    auto x = max_x(i, j);
    return {i, j, value, x};
  }

  friend std::ostream &operator<<(std::ostream &os, const Matrix_Diff &obj)
  {
    const int height = obj.max_diff.Height();
    const int width = obj.max_diff.Width();
    ASSERT_EQUAL(height, obj.max_x.Height());
    ASSERT_EQUAL(width, obj.max_x.Width());

    const auto label = "";
    os << El::DimsString(obj.max_diff, label);
    os << "\n  [\n  ";
    for(int i = 0; i < height; ++i)
      {
        for(int j = 0; j < width; ++j)
          {
            // os << " " << matrix.Get(i, j);
            os << " (max=" << obj.max_diff(i, j) << ", x=" << obj.max_x(i, j)
               << ")";
          }
        os << "\n  ";
      }
    os << "]";

    auto [i, j, max_diff, x] = obj.get_max_element_and_x();
    os << "\n  max_diff=" << max_diff << " at (" << i << "," << j
       << ") x=" << x;
    return os;
  }
};

// For each poly_vec element of PMP matrix,
// compute product poly_vec.z and evaluate at sample points (including prefactor)
// => result is vector of matrices, length = num_points
[[nodiscard]] std::vector<El::Matrix<El::BigFloat>>
sample_pmp_matrix(const Polynomial_Vector_Matrix &matrix,
                  const El::Matrix<El::BigFloat> &z)
{
  std::vector<El::Matrix<El::BigFloat>> result;
  const size_t num_points = matrix.sample_points.size();
  for(int point_index = 0; point_index < num_points; ++point_index)
    {
      auto &matrix_sample = result.emplace_back(matrix.polynomials.Height(),
                                                matrix.polynomials.Width());
      El::Zero(matrix_sample);
      const auto &x = matrix.sample_points.at(point_index);
      const auto &scale = matrix.sample_scalings.at(point_index);
      for(int i = 0; i < matrix.polynomials.Height(); ++i)
        for(int j = 0; j < matrix.polynomials.Width(); ++j)
          {
            const auto &poly_vec = matrix.polynomials.CRef(i, j);

            ASSERT_EQUAL(poly_vec.size(), z.Height());
            El::BigFloat &value = matrix_sample.Ref(i, j);
            for(int z_index = 0; z_index < z.Height(); ++z_index)
              {
                const auto &poly = poly_vec.at(z_index);
                value += z(z_index, 0) * scale * poly(x);
              }
          }
    }
  return result;
}

// Attempting to measure difference of two sample matrices
Matrix_Diff matrix_diff(const std::vector<El::Matrix<El::BigFloat>> &left,
                        const std::vector<El::Matrix<El::BigFloat>> &right,
                        const std::vector<El::BigFloat> &sample_points)
{
  size_t num_points = sample_points.size();
  ASSERT_EQUAL(left.size(), num_points);
  ASSERT_EQUAL(right.size(), num_points);

  const auto height = left.at(0).Height();
  const auto width = left.at(0).Width();

  Matrix_Diff result;
  result.resize(height, width);
  for(int i = 0; i < height; ++i)
    for(int j = 0; j < width; ++j)
      {
        auto &max_diff = result.max_diff.Ref(i, j);
        auto &max_x = result.max_x.Ref(i, j);

        for(size_t index = 0; index < num_points; ++index)
          {
            auto curr_diff
              = El::Abs(left.at(index)(i, j) - right.at(index)(i, j));
            if(curr_diff > max_diff)
              {
                max_diff = curr_diff;
                max_x = sample_points.at(index);
              }
          }
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
std::vector<std::tuple<fs::path, fs::path, Matrix_Diff>>
diff_pmp(const Polynomial_Matrix_Program &left_pmp,
         const Polynomial_Matrix_Program &right_pmp,
         const El::Matrix<El::BigFloat> &z)
{
  ASSERT_EQUAL(left_pmp.matrices.size(), right_pmp.matrices.size());
  ASSERT(left_pmp.matrix_index_local_to_global
         == right_pmp.matrix_index_local_to_global);

  const size_t num_matrices = left_pmp.num_matrices;
  // block_index, left, rigt, result
  std::vector<std::tuple<fs::path, fs::path, Matrix_Diff>> result(
    num_matrices);

  for(size_t local_index = 0; local_index < left_pmp.matrices.size();
      ++local_index)
    {
      const auto &left_pvm = left_pmp.matrices.at(local_index);
      const auto &right_pvm = right_pmp.matrices.at(local_index);

      const auto &left_path = left_pmp.block_paths.at(local_index);
      const auto &right_path = right_pmp.block_paths.at(local_index);

      result.at(left_pmp.matrix_index_local_to_global.at(local_index))
        = {left_path, right_path, diff_pvm(left_pvm, right_pvm, z)};
    }

  const int rank = El::mpi::Rank();
  std::vector<int> block_index_to_rank(num_matrices, El::mpi::Size());
  for(auto &block_index : left_pmp.matrix_index_local_to_global)
    {
      block_index_to_rank.at(block_index) = rank;
    }
  const auto comm = El::mpi::COMM_WORLD;
  El::mpi::AllReduce(block_index_to_rank.data(), block_index_to_rank.size(),
                     El::mpi::MIN, comm);

  // Send all results to rank=0 and print
  for(size_t block_index = 0; block_index < num_matrices; ++block_index)
    {
      const int from = block_index_to_rank.at(block_index);
      const int to = 0;

      if(rank != to && rank != from)
        continue;

      auto &[left_path, right_path, matrix_diff] = result.at(block_index);
      copy_matrix(matrix_diff.max_diff, from, to, comm);
      copy_matrix(matrix_diff.max_x, from, to, comm);

      // Copy paths to rank=0
      for(auto *path : {&left_path, &right_path})
        {
          auto path_string = path->string();
          if(rank == from)
            {
              const size_t path_size = path_string.size();
              El::mpi::Send<size_t>(path_size, to, El::mpi::COMM_WORLD);
              MPI_Send(path_string.data(), path_size, MPI_CHAR, to, 0,
                       MPI_COMM_WORLD);
            }
          if(rank == to)
            {
              const auto path_size
                = El::mpi::Recv<size_t>(from, El::mpi::COMM_WORLD);
              path_string.resize(path_size);
              MPI_Recv(path_string.data(), path_size, MPI_CHAR, from, 0,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              *path = path_string;
            }
        }
    }
  return result;
}

void compare_pmp(const Polynomial_Matrix_Program &left_pmp,
                 const Polynomial_Matrix_Program &right_pmp,
                 const El::Matrix<El::BigFloat> &z,
                 const fs::path &output_path, const bool sort_stdout,
                 const Verbosity verbosity)
{
  auto diff = diff_pmp(left_pmp, right_pmp, z);
  //Print results
  if(El::mpi::Rank() != 0)
    return;

  std::vector<std::tuple<size_t, fs::path, fs::path, Matrix_Diff>>
    sorted_result;
  {
    for(size_t block_index = 0; block_index < diff.size(); ++block_index)
      {
        const auto &[left_path, right_path, value] = diff.at(block_index);
        sorted_result.emplace_back(block_index, left_path, right_path, value);
      }
    if(sort_stdout)
      {
        auto get_key =
          [](const std::tuple<size_t, fs::path, fs::path, Matrix_Diff> &item) {
            const auto &[block_index, left, right, value] = item;
            const auto &[i, j, max_value, x] = value.get_max_element_and_x();
            return std::make_tuple(max_value, block_index);
          };
        std::sort(sorted_result.begin(), sorted_result.end(),
                  [&get_key](const auto &a, const auto &b) {
                    return get_key(a) < get_key(b);
                  });
      }
  }

  // Print to stdout
  if(verbosity >= Verbosity::regular)
    {
      El::Output("======================================");
      El::Output("Comparison results",
                 sort_stdout ? " sorted by max diff:" : ":");
      El::Output("======================================");
      for(const auto &[block_index, left_path, right_path, value] :
          sorted_result)
        {
          El::Output("block_index=", block_index);
          El::Output("  left=", left_path);
          El::Output("  right=", right_path);
          El::Output("  diff=", value);
          El::Output("--------------------------------------");
        }
      El::Output("======================================");
    }

  // write to JSON
  if(!output_path.empty())
    {
      if(verbosity >= Verbosity::regular)
        El::Output("Writing JSON output to ", weakly_canonical(output_path).string());

      if(output_path.has_parent_path())
        fs::create_directories(output_path.parent_path());
      std::ofstream os(output_path);
      ASSERT(os.good(), "Cannot open --output=", output_path);
      rapidjson::OStreamWrapper osw(os);
      rapidjson::Writer writer(osw);

      writer.StartArray();
      for(size_t block_index = 0; block_index < diff.size(); ++block_index)
        {
          writer.StartObject();
          const auto &[left_path, right_path, value] = diff.at(block_index);
          const auto &[max_i, max_j, max_value, max_x]
            = value.get_max_element_and_x();
          writer.Key("block_index");
          writer.Uint(block_index);
          writer.Key("left_path");
          writer.String(left_path.c_str());
          writer.Key("right_path");
          writer.String(right_path.c_str());
          writer.Key("max");
          writer.Double(static_cast<double>(max_value));
          writer.Key("i");
          writer.Int(max_i);
          writer.Key("j");
          writer.Int(max_j);
          writer.Key("x");
          writer.Double(static_cast<double>(max_x));

          writer.Key("diff_matrix");
          writer.StartArray();
          for(int i = 0; i < value.max_diff.Height(); ++i)
            {
              writer.StartArray();
              for(int j = 0; j < value.max_diff.Width(); ++j)
                {
                  writer.Double(static_cast<double>(value.max_diff(i, j)));
                }
              writer.EndArray();
            }
          writer.EndArray();

          writer.Key("x_matrix");
          writer.StartArray();
          for(int i = 0; i < value.max_x.Height(); ++i)
            {
              writer.StartArray();
              for(int j = 0; j < value.max_x.Width(); ++j)
                {
                  writer.Double(static_cast<double>(value.max_x(i, j)));
                }
              writer.EndArray();
            }
          writer.EndArray();

          writer.EndObject();
        }
      writer.EndArray();
    }
}

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      int precision;
      fs::path left_pmp_path, right_pmp_path;
      fs::path z_path;
      fs::path output_path;
      bool sort_stdout = false;
      Verbosity verbosity = Verbosity::debug;

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
        options.add_options()("output,o", po::value<fs::path>(&output_path),
                              "Directory to place output");
        options.add_options()("sort,s", po::bool_switch(&sort_stdout),
                              "Sort (by max_diff) results printed to stdout.");
        options.add_options()(
          "verbosity", po::value<Verbosity>(&verbosity),
          "Verbosity (0,1,2,3 or none,regular,debug,trace)");

        po::positional_options_description positional;
        positional.add("precision", 1);
        positional.add("left", 1);
        positional.add("right", 1);
        positional.add("z", 1);
        positional.add("output", 1);

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
      compare_pmp(left_pmp, right_pmp, z, output_path, sort_stdout, verbosity);
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
