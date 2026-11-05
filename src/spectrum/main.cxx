#include "Zeros.hxx"
#include "pmp/PMP_Info.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "sdpb_util/Boost_Float.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void handle_arguments(const int &argc, char **argv,
                      std::optional<Boost_Float> &threshold,
                      El::BigFloat &max_zero, El::BigFloat &min_zero_distance,
                      fs::path &pmp_info_path, fs::path &solution_dir,
                      fs::path &c_minus_By_path, fs::path &output_path,
                      bool &need_lambda,
                      std::optional<El::BigFloat> &min_eigenvalue_ratio,
                      Verbosity &verbosity);

void set_default_parameters(const std::filesystem::path &solution_dir,
                            const PMP_Info &pmp_info, const bool &need_lambda,
                            const Verbosity &verbosity, Timers &timers,
                            std::optional<Boost_Float> &threshold,
                            std::optional<El::BigFloat> &min_eigenvalue_ratio);

PMP_Info
read_pmp_info(const std::filesystem::path &input_path, Timers &timers);

std::vector<El::Matrix<El::BigFloat>>
read_c_minus_By(const std::filesystem::path &input_path,
                const PMP_Info &pmp_info, Timers &timers);

std::vector<El::Matrix<El::BigFloat>>
read_x(const fs::path &solution_path, const PMP_Info &pmp_info,
       Timers &timers);

std::vector<Zeros>
compute_spectrum(const PMP_Info &pmp_info,
                 const std::vector<El::Matrix<El::BigFloat>> &c_minus_By,
                 const std::optional<std::vector<El::Matrix<El::BigFloat>>> &x,
                 const Boost_Float &threshold, const El::BigFloat &max_zero,
                 const El::BigFloat &min_zero_distance,
                 const bool &need_lambda,
                 const std::optional<El::BigFloat> &min_eigenvalue_ratio,
                 const Verbosity &verbosity,
                 const std::filesystem::path &spectrum_output_path,
                 Timers &timers);

void write_spectrum(const fs::path &output_path,
                    const std::vector<Zeros> &zeros_blocks,
                    const PMP_Info &pmp_info, Timers &timers);

void create_profiling_dir(const fs::path &spectrum_output_path);
void write_profiling(const fs::path &spectrum_output_path, Timers &timers);

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      std::optional<Boost_Float> threshold;
      El::BigFloat max_zero;
      El::BigFloat min_zero_distance;
      fs::path pmp_info_path, solution_dir, output_path, c_minus_By_path;
      bool need_lambda;
      std::optional<El::BigFloat> min_eigenvalue_ratio;
      Verbosity verbosity;
      handle_arguments(argc, argv, threshold, max_zero, min_zero_distance,
                       pmp_info_path, solution_dir, c_minus_By_path,
                       output_path, need_lambda, min_eigenvalue_ratio,
                       verbosity);

      // Print command line
      if(verbosity >= Verbosity::debug && El::mpi::Rank() == 0)
        {
          std::vector<std::string> arg_list(argv, argv + argc);
          for(const auto &arg : arg_list)
            std::cout << arg << " ";
          std::cout << std::endl;
        }

      Timers timers(env, verbosity);
      Scoped_Timer timer(timers, "spectrum");
      const auto pmp_info = read_pmp_info(pmp_info_path, timers);

      // Set --threshold and --minEigenvalueRatio to sqrt(dualityGap), if needed.
      set_default_parameters(solution_dir, pmp_info, need_lambda, verbosity,
                             timers, threshold, min_eigenvalue_ratio);

      std::optional<std::vector<El::Matrix<El::BigFloat>>> x;
      if(need_lambda)
        x.emplace(read_x(solution_dir, pmp_info, timers));

      const auto c_minus_By
        = read_c_minus_By(c_minus_By_path, pmp_info, timers);

      // Create directory spectrum.json.profiling/
      if(verbosity >= Verbosity::debug)
        create_profiling_dir(output_path);

      ASSERT(threshold.has_value(), "--threshold must be specified!");
      const auto zeros_blocks = compute_spectrum(
        pmp_info, c_minus_By, x, threshold.value(), max_zero,
        min_zero_distance, need_lambda, min_eigenvalue_ratio, verbosity,
        output_path, timers);

      write_spectrum(output_path, zeros_blocks, pmp_info, timers);

      // Write profiling data
      if(verbosity >= Verbosity::debug)
        write_profiling(output_path, timers);
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
