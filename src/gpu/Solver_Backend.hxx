// GPU backend for SDPB: the boundary between SDPB's driver and the arithmetic
// (design sketch; not yet compiled into sdpb).
//
// The driver (SDP_Solver::run: SIGTERM, checkpoint timing, termination and
// jump tests, iterations.json, c_minus_By, outputs) stays SDPB's code and is
// shared by both backends.  Each iteration it calls the four operations below
// and sees only scalars, as El::BigFloat (exact conversions from device floats).
// The CPU backend wraps today's functions unchanged; the GPU backend runs them
// chunk-major on the device, as SDPviBe's solver.py does (Antoine Tilloy,
// https://github.com/AntoineTilloy/SDPviBe, the prototype of this backend).
#pragma once

#include <El.hpp>

#include <filesystem>
#include <string>

class SDP;
class Block_Info;
class Solver_Parameters;
class SDP_Solver;

namespace sdpb_gpu
{
  // What compute_objectives, compute_*_residues_and_error, and the start of
  // step() produce: everything compute_feasible_and_termination and
  // print_iteration read.
  struct Residue_Scalars
  {
    El::BigFloat primal_objective, dual_objective, duality_gap;
    El::BigFloat primal_error_P, primal_error_p, dual_error;
    El::BigFloat mu;                      // Tr(XY) / total_psd_rows
    El::BigFloat max_block_cond_number;   // X, Y Cholesky factors
    std::string max_block_cond_number_name;
  };

  struct Factor_Scalars
  {
    El::BigFloat Q_cond_number, max_block_cond_number; // Schur blocks, Q
    std::string max_block_cond_number_name;
  };

  struct Direction_Scalars
  {
    El::BigFloat R_error;           // predictor only
    El::BigFloat frobenius_product; // Tr((X + dX)(Y + dY)), predictor only
  };

  class Solver_Backend
  {
  public:
    virtual ~Solver_Backend() = default;

    // the SDP (parsed by SDPB's own readers) and the starting point (SDPB's
    // initial matrices or a loaded checkpoint), converted once to the device
    virtual void load(const SDP &sdp, const Block_Info &block_info,
                      const SDP_Solver &initial)
      = 0;

    // pending update (x, X, y, Y) += (alpha_P dx, alpha_P dX, alpha_D dy,
    // alpha_D dY) from the last step; Cholesky factors of X, Y; bilinear
    // pairings; dual and primal residues; objectives; errors; mu
    virtual Residue_Scalars residues(const Solver_Parameters &p) = 0;

    // Schur blocks S_j, their factors, P = L^-1 B (kept when it fits), Q and
    // its factor (Gram_Accumulator, then a tiled Cholesky)
    virtual Factor_Scalars factor() = 0;

    // predictor (corrector = false) or corrector direction for
    // centering parameter beta
    virtual Direction_Scalars
    direction(const El::BigFloat &beta, const El::BigFloat &mu, bool corrector)
      = 0;

    // primal and dual step lengths (stepLengthReduction gamma); the update is
    // applied lazily by the next residues() (SDPviBe's order of passes)
    virtual std::pair<El::BigFloat, El::BigFloat>
    step_lengths(const El::BigFloat &gamma) = 0;
    virtual void set_steps(const El::BigFloat &alpha_P,
                           const El::BigFloat &alpha_D)
      = 0;

    // x, X, y, Y back to SDPB's types (checkpoints, outputs, c_minus_By)
    virtual void get_iterate(SDP_Solver &out) = 0;
  };
}
