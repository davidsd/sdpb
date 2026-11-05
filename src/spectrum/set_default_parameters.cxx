#include "pmp/PMP_Info.hxx"
#include "sdp_solve/read_text_block.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <boost/algorithm/string.hpp>

namespace fs = std::filesystem;

namespace
{
  El::BigFloat read_duality_gap_impl(const fs::path &solution_dir)
  {
    const auto out_txt_path = solution_dir / "out.txt";
    std::ifstream is(out_txt_path);
    ASSERT(is.good(), "Failed to open ", out_txt_path);
    std::string line;
    while(std::getline(is, line))
      {
        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of("=;"));
        for(auto &token : tokens)
          boost::trim(token);

        if(tokens.empty())
          break;
        if(tokens[0] != "dualityGap")
          continue;
        const auto duality_gap = El::BigFloat(tokens[1]);
        return duality_gap;
      }
    RUNTIME_ERROR("Cannot find dualityGap in ", out_txt_path);
  }

  El::BigFloat read_duality_gap(const fs::path &solution_dir, Timers &timers)
  {
    Scoped_Timer timer(timers, "read_duality_gap");
    El::BigFloat duality_gap;
    if(El::mpi::Rank() == 0)
      duality_gap = read_duality_gap_impl(solution_dir);
    El::mpi::Broadcast(duality_gap, 0, El::mpi::COMM_WORLD);
    return duality_gap;
  }
}

void set_default_parameters(
  const std::filesystem::path &solution_dir, const PMP_Info &pmp_info,
  const bool &need_lambda, const Verbosity &verbosity, Timers &timers,
  std::optional<Boost_Float> &threshold,
  std::optional<El::BigFloat> &min_eigenvalue_ratio)
{
  const bool need_threshold = !threshold.has_value();

  const bool need_min_eigenvalue_ratio = [&] {
    if(min_eigenvalue_ratio.has_value())
      return false;
    // --minEigenvalueRatio needed only in compute_lambda
    if(!need_lambda)
      return false;
    // --minEigenvalueRatio is not needed for 1x1 blocks.
    size_t max_block_dim = 0;
    for(const auto &block : pmp_info.blocks)
      max_block_dim = std::max(max_block_dim, block.dim);
    max_block_dim
      = El::mpi::AllReduce(max_block_dim, El::mpi::MAX, El::mpi::COMM_WORLD);
    return max_block_dim > 1;
  }();

  if(!need_threshold && !need_min_eigenvalue_ratio)
    return;

  // Read dualityGap
  const auto sqrt_duality_gap
    = El::Sqrt(read_duality_gap(solution_dir, timers));
  if(El::mpi::Rank() == 0)
    {
      ASSERT(sqrt_duality_gap >= 0 && sqrt_duality_gap <= 1,
             "sqrt(dualityGap)=", sqrt_duality_gap,
             " is expected to be in range [0,1]");
    }
  // Wait for checks performed at rank=0
  El::mpi::Barrier();

  // --threshold = sqrt(dualityGap)
  if(need_threshold)
    {
      threshold.emplace(to_Boost_Float(sqrt_duality_gap));
      if(El::mpi::Rank() == 0 && verbosity >= Verbosity::regular)
        {
          El::Output("Setting --threshold to sqrt(dualityGap)=",
                     sqrt_duality_gap);
        }
    }

  // --minEigenvalueRatio=sqrt(dualityGap)
  if(need_min_eigenvalue_ratio)
    {
      min_eigenvalue_ratio.emplace(sqrt_duality_gap);
      if(El::mpi::Rank() == 0 && verbosity >= Verbosity::regular)
        {
          El::Output("Setting --minEigenvalueRatio to sqrt(dualityGap)=",
                     sqrt_duality_gap);
        }
    }
}
