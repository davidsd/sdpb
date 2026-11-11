#pragma once

// sdpb --gpu on several MPI ranks (one GPU each): every rank owns whole
// GPU chunks of blocks, in place of SDPB's block-to-group mapping.
// Compiled only with --enable-gpu.

#include "SDPB_Parameters.hxx"
#include "gpu/solver/Gpu_Layout.hxx"
#include "sdp_solve/Block_Info.hxx"
#include "sdpa_solve/Block_Info.hxx"
#include "sdpb_util/Environment.hxx"
#include "sdpb_util/assert.hxx"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

extern char **environ;

// The SDPB_GPU_* settings are read from each process's environment, so a
// launcher that doesn't pass them to every rank would give each rank its own
// placement and arithmetic. Refuse that.
inline void check_gpu_settings_agree()
{
  std::vector<std::string> settings;
  for(char **e = environ; *e != nullptr; ++e)
    if(std::strncmp(*e, "SDPB_GPU_", 9) == 0)
      settings.emplace_back(*e);
  std::sort(settings.begin(), settings.end());
  uint64_t hash = 1469598103934665603ULL; // FNV-1a
  for(const auto &setting : settings)
    for(const unsigned char c : setting + '\n')
      hash = (hash ^ c) * 1099511628211ULL;
  uint64_t lo = 0, hi = 0;
  MPI_Allreduce(&hash, &lo, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&hash, &hi, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  ASSERT(lo == hi, "--gpu: the SDPB_GPU_* environment settings differ "
                   "between MPI ranks; pass the same ones to every rank "
                   "(e.g. mpirun -x).");
}

// The block placement is written next to the checkpoints (gpu_layout), and a
// checkpoint is loaded only with the placement it was written with: the
// binary checkpoint holds each rank's local blocks in order, so another
// placement would misread it (or, if the shapes happen to match, read other
// blocks' data).
inline void check_and_record_gpu_layout(const SDPB_Parameters &parameters,
                                        const std::vector<int> &block_owner)
{
  namespace fs = std::filesystem;
  std::ostringstream layout;
  layout << "# sdpb --gpu block placement: block rank\n";
  for(size_t block = 0; block < block_owner.size(); ++block)
    layout << block << " " << block_owner[block] << "\n";
  const auto &in = parameters.solver.checkpoint_in;
  int verdict = 0; // 1: no gpu_layout with the checkpoint, 2: another one
  if(El::mpi::Rank() == 0)
    {
      if(fs::exists(in / "checkpoint.json") || fs::exists(in / "checkpoint.0"))
        {
          std::ifstream stream(in / "gpu_layout");
          std::stringstream recorded;
          recorded << stream.rdbuf();
          verdict = !stream.is_open()                 ? 1
                    : recorded.str() != layout.str() ? 2
                                                     : 0;
        }
      if(verdict == 0)
        {
          const auto &out = parameters.solver.checkpoint_out;
          fs::create_directories(out);
          std::ofstream(out / "gpu_layout") << layout.str();
        }
    }
  El::mpi::Broadcast(verdict, 0, El::mpi::COMM_WORLD);
  ASSERT(verdict != 1, "--gpu on ", El::mpi::Size(), " ranks: the checkpoint in ",
         in, " has no gpu_layout, so it may come from another block placement "
         "(a CPU run, other SDPB_GPU_CHUNK_* settings or rank counts).");
  ASSERT(verdict != 2, "--gpu on ", El::mpi::Size(), " ranks: the checkpoint in ",
         in, " was written with another GPU block placement (its gpu_layout "
         "differs: other SDPB_GPU_CHUNK_* settings or rank counts).");
}

inline void
assign_blocks_to_gpus(const SDPB_Parameters &parameters, Block_Info &block_info)
{
  check_gpu_settings_agree();
  ASSERT(block_info.num_nodes() == 1,
         "The distributed GPU prototype currently requires one node");
  std::vector<sdpb_gpu::Block_Shape> shapes(block_info.dimensions.size());
  for(size_t j = 0; j < shapes.size(); ++j)
    shapes[j] = {int(block_info.dimensions[j]), int(block_info.num_points[j]),
                 int(block_info.get_bilinear_bases_height(j, 0)),
                 int(block_info.get_bilinear_bases_height(j, 1))};
  const char *cap = std::getenv("SDPB_GPU_CHUNK_BLOCKS");
  const auto chunks = sdpb_gpu::gpu_chunks(shapes, cap ? std::atoi(cap) : 16);
  const auto owners = sdpb_gpu::gpu_chunk_owners(chunks, El::mpi::Size());
  std::vector<int> block_owner(shapes.size(), -1);
  for(size_t c = 0; c < chunks.size(); ++c)
    for(const int j : chunks[c].blocks)
      block_owner[j] = owners[c];
  block_info.allocate_single_rank_blocks(block_owner);
  check_and_record_gpu_layout(parameters, block_owner);
}

inline void
assign_blocks_to_gpus(const SDPB_Parameters &, Sdpb::Sdpa::Block_Info &)
{
  RUNTIME_ERROR("--gpu is implemented for polynomial matrix programs only, "
                "not for SDPA input.");
}
