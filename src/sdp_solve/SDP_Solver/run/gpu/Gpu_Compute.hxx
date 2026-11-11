// GPU backend for SDPB: what SDP_Solver::run calls in place of SDPB's own
// per-iteration functions when --gpu is given.  Converts between El::BigFloat and
// device floats; the arithmetic is in src/gpu/solver/Gpu_Solver (based on SDPviBe
// by Antoine Tilloy, https://github.com/AntoineTilloy/SDPviBe, the prototype of
// this backend; docs/GPU_backend_credits.md).
#pragma once

#include "sdp_solve/Block_Info.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdp_solve/Solver_Parameters.hxx"

#include <El.hpp>

#include <filesystem>
#include <memory>
#include <string>

struct SDP_Solver;

class Gpu_Compute
{
public:
  Gpu_Compute(const Solver_Parameters &parameters, const Block_Info &block_info,
              const SDP &sdp, const SDP_Solver &solver);
  ~Gpu_Compute();

  // Limbs L of the device floats for --precision (32 L - 2 mantissa bits).
  // Throws, naming --precision, if this build has no kernels for it, so that
  // sdpb can refuse it at start-up.
  static int limbs_for_precision(size_t precision);

  // Objectives, errors: as compute_objectives and compute_*_residues_and_error,
  // into solver's scalars
  void residues(const SDP &sdp, SDP_Solver &solver);

  // As SDP_Solver::step
  void
  step(const Solver_Parameters &parameters, size_t total_psd_rows,
       bool is_primal_and_dual_feasible, SDP_Solver &solver, El::BigFloat &mu,
       El::BigFloat &beta_corrector, El::BigFloat &primal_step_length,
       El::BigFloat &dual_step_length, bool &terminate_now,
       El::BigFloat &Q_cond_number, El::BigFloat &max_block_cond_number,
       std::string &max_block_cond_number_name);

  // x, X, y, Y back into solver (checkpoints, outputs, c_minus_By)
  void get_iterate(SDP_Solver &solver);

  void write_profile(const std::filesystem::path &path) const;

private:
  struct Impl;
  std::unique_ptr<Impl> p_;
};
