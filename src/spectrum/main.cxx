#include "Zeros.hxx"
#include "pmp/PMP_Info.hxx"
#include "sdp_solve/sdp_solve.hxx"

#include <filesystem>

namespace fs = std::filesystem;

void handle_arguments(const int &argc, char **argv, El::BigFloat &threshold,
                      fs::path &pmp_info_path, fs::path &solution_dir,
                      fs::path &c_minus_By_path, fs::path &output_path,
                      bool &need_lambda, Verbosity &verbosity);

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
                 const El::BigFloat &threshold, const bool &need_lambda,
                 const Verbosity &verbosity, Timers &timers);

void write_spectrum(const fs::path &output_path,
                    const std::vector<Zeros> &zeros_blocks,
                    const PMP_Info &pmp_info, Timers &timers);

void write_profiling(const fs::path &spectrum_output_path, Timers &timers);

int main(int argc, char **argv)
{
  Environment env(argc, argv);

  try
    {
      El::BigFloat threshold;
      fs::path pmp_info_path, solution_dir, output_path, c_minus_By_path;
      bool need_lambda;
      Verbosity verbosity;
      handle_arguments(argc, argv, threshold, pmp_info_path, solution_dir,
                       c_minus_By_path, output_path, need_lambda, verbosity);

      // Print command line
      if(verbosity >= Verbosity::debug && El::mpi::Rank() == 0)
        {
          std::vector<std::string> arg_list(argv, argv + argc);
          for(const auto &arg : arg_list)
            std::cout << arg << " ";
          std::cout << std::endl;
        }

      // TODO use timers, print profiling data for --verbosity=debug
      Timers timers(env, verbosity);
      Scoped_Timer timer(timers, "spectrum");
      const auto pmp_info = read_pmp_info(pmp_info_path, timers);

      std::optional<std::vector<El::Matrix<El::BigFloat>>> x;
      if(need_lambda)
        x.emplace(read_x(solution_dir, pmp_info, timers));

      const auto c_minus_By
        = read_c_minus_By(c_minus_By_path, pmp_info, timers);
      const auto zeros_blocks = compute_spectrum(
        pmp_info, c_minus_By, x, threshold, need_lambda, verbosity, timers);

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
