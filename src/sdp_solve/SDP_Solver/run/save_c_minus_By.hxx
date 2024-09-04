#pragma once

#include "sdp_solve/Block_Vector.hxx"
#include "sdp_solve/SDP.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>

#include <filesystem>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/rapidjson.h>
#include <rapidjson/writer.h>

// Save vector c - B.y
// TODO split calculation and writing
// TODO reuse calculation in compute_dual_residues_and_error()
inline void save_c_minus_By(const std::filesystem::path &path,
                            const Block_Info &block_info, const SDP &sdp,
                            const Block_Vector &y, const Verbosity &verbosity,
                            Timers &timers)
{
  Scoped_Timer timer(timers, "save_c_minus_By");

  // dist_c_minus_By := alpha * B.y + beta * c == c - B.y
  Block_Vector dist_c_minus_By = sdp.primal_objective_c;
  {
    const El::BigFloat alpha(-1);
    const El::BigFloat beta(1);
    for(int index = 0; index < dist_c_minus_By.blocks.size(); ++index)
      {
        const auto &B_block = sdp.free_var_matrix.blocks.at(index);
        El::Gemm(El::Orientation::NORMAL, El::Orientation::NORMAL, alpha,
                 B_block, y.blocks.at(index), beta,
                 dist_c_minus_By.blocks.at(index));
      }
  }

  // Copy from distributed BlockVector blocks to local vector<Matrix> on rank=0
  // TODO extract to a separate function
  const int rank = El::mpi::Rank();
  const size_t num_blocks = block_info.dimensions.size();
  // To be filled on rank=0
  std::vector<El::Matrix<El::BigFloat>> c_minus_By(num_blocks);

  for(size_t block_index = 0; block_index < num_blocks; ++block_index)
    {
      // Each block is copied to the root of its MPI group and stored in send_matrix
      std::optional<El::Matrix<El::BigFloat>> send_matrix;

      // ranks for Send/Recv:
      int from = -1;
      const int to = 0;

      for(size_t local_block_index = 0;
          local_block_index < block_info.block_indices.size();
          ++local_block_index)
        {
          if(block_info.block_indices.at(local_block_index) != block_index)
            continue;
          // Copy all elements to the root of the current MPI Group
          const auto &block = dist_c_minus_By.blocks.at(local_block_index);
          El::DistMatrix<El::BigFloat, El::CIRC, El::CIRC> block_circ_circ(
            block);

          // NB: Root() and Grid().Rank() return ranks inside MPI group,
          // whereas 'from' is a global rank
          if(block_circ_circ.Grid().Rank() == block_circ_circ.Root())
            {
              from = rank;
              send_matrix.emplace(block_circ_circ.LockedMatrix());
            }
        }

      from = El::mpi::AllReduce(from, El::mpi::MAX, El::mpi::COMM_WORLD);
      if(rank == to && rank == from)
        {
          // No need for MPI in this case
          c_minus_By.at(block_index) = send_matrix.value();
          continue;
        }

      if(rank == from)
        {
          ASSERT_EQUAL(send_matrix->Width(), 1, DEBUG_STRING(from),
                       DEBUG_STRING(block_index),
                       El::DimsString(*send_matrix, " send_matrix"));
          El::mpi::Send(send_matrix->Height(), to, El::mpi::COMM_WORLD);
          El::Send(send_matrix.value(), El::mpi::COMM_WORLD, to);
        }
      if(rank == to)
        {
          auto &block = c_minus_By.at(block_index);
          int height = El::mpi::Recv<int>(from, El::mpi::COMM_WORLD);
          block.Resize(height, 1);
          El::Recv(block, El::mpi::COMM_WORLD, from);
        }
    }

  // Write block vector to file.
  // Example:
  // If there are two blocks, each having three element,
  // then the JSON output will look like:
  // {"c_minus_By" : [["111","222","333], ["44","55","66"]]}
  {
    if(rank != 0)
      return;

    if(verbosity >= Verbosity::debug)
      {
        El::Output("Saving (c-B.y) to: ", path);
      }

    if(path.has_parent_path())
      create_directories(path.parent_path());

    std::ofstream os(path);
    if(!os.good())
      {
        PRINT_WARNING("Cannot write (c-B.y) to ", path);
        return;
      }
    rapidjson::OStreamWrapper osw(os);
    rapidjson::Writer writer(osw);
    writer.StartObject();
    writer.Key("c_minus_By");
    writer.StartArray();

    // reusable stream for BigFloats
    std::stringstream ss;
    set_stream_precision(ss);
    for(const auto &block : c_minus_By)
      {
        ASSERT_EQUAL(block.Width(), 1);
        writer.StartArray();
        for(int i = 0; i < block.Height(); ++i)
          {
            ss.str({});
            ss << block.Get(i, 0);
            writer.String(ss.str().c_str());
          }
        writer.EndArray();
      }
    writer.EndArray();
    writer.EndObject();
    ASSERT(writer.IsComplete());
  }
}