#pragma once

#include "sdpb/SDPB_Parameters.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdpb_util/Environment.hxx"

#include <gmp.h>

#include <string>
#include <vector>

// Resident mode, sdpb --controlDir <dir>.
//
// Load the SDP once, then execute commands from the control directory
// until 'quit': solve from the current iterate, replace blocks or objectives
// in memory, set the coefficients of a parametric block, save and restore
// iterates, write checkpoints. See Readme.md in this directory.
//
// startup_args: the command line of the resident process (without argv[0]).
// Options missing from a solve command are taken from it.
// startup_precision: GMP's precision when the command line was parsed. Solve
// commands are parsed at it too, so that their BigFloat options (thresholds,
// centering parameters, ...) have the same values as in an sdpb run.
void run_resident(const Environment &env, const SDPB_Parameters &parameters,
                  const std::vector<std::string> &startup_args,
                  mp_bitcnt_t startup_precision, const Block_Info &block_info);
