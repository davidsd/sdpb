#include "../BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdp_solve/SDP_Solver/run/bigint_syrk/fmpz/Fmpz_BigInt.hxx"
#include "restore_bigint_from_residues.hxx"
#include "sdpb_util/block_mapping/MPI_Comm_Wrapper.hxx"

// Restore output matrix from residues (output_residues_window) and synchronize it for all nodes.
//
// Initially, each node n has its own matrix Q_n, stored as residues in output_residues_window.
// The function restores Q_n from residues, sums Q_n over all nodes
// and writes the resulting Q into a global DistMatrix output.
//
// Algorithm:
// 1. The rank that owns the element Q[i,j] sets Q[i,j] = Q_n[i,j] restored from residues on its node.
// 2. Accumulate Q_n[i,j] from other nodes to Q[i,j] according to the following scheme:
//   for offset = 1..num_nodes-1:
//   - Each node (n) restores all Q_n[i,j] for [i,j] owned by node (n+offset) from residues,
//     and sends them to node (n+offset). Implemented via MPI_Sendrecv.
//   - Each node updates its own Q[i,j] with data received from node (n-offset).
//
// All ranks on a node have access to each Q_n[i,j], so how do we decide
//   which rank should send it?
// If rank r owns the given element Q[i,j], then it will receive contributions
//   from other nodes from the processes having the same rank within a node.
// For example, for 3 nodes with 128 cores,
//   rank 0 will communicate only with ranks 128 and 256,
//   rank 1 - with ranks 129 and 257, and so on.
//   Communicator for such group, e.g. (0,128,256), (1,129,257) etc.,
//     is called reduce_comm in the code.
//
// All [i,j] owned by one rank are combined into a single MPI message.
// Thus, each rank has send and receive buffers of size ~ #(Q) / num_ranks.
// Each rank performs (num_nodes - 1) send/recv operations.
//
void BigInt_Shared_Memory_Syrk_Context::restore_and_reduce(
  std::optional<El::UpperOrLower> uplo, El::DistMatrix<El::BigFloat> &output,
  Timers &timers)
{
  Scoped_Timer timer(timers, "restore_and_reduce");

  auto skip_element = [&uplo](auto row, auto column) {
    if(uplo.has_value())
      {
        if(*uplo == El::UPPER && column < row)
          return true;
        if(*uplo == El::LOWER && row < column)
          return true;
      }
    return false;
  };

  const auto output_comm = output.DistComm();
  const auto height = output.Height();
  const auto width = output.Width();
  ASSERT(height <= output_residues_window->height, DEBUG_STRING(height),
         DEBUG_STRING(output_residues_window->height));
  ASSERT(width <= output_residues_window->width, DEBUG_STRING(width),
         DEBUG_STRING(output_residues_window->width));

  Fmpz_BigInt bigint_value;
  std::vector<mp_limb_t> residues_buffer_temp(comb.num_primes);

  {
    Scoped_Timer local_restore_timer(timers, "local_restore");
    for(int i = 0; i < height; ++i)
      for(int j = 0; j < width; ++j)
        {
          if(skip_element(i, j))
            continue;
          if(!output.IsLocal(i, j))
            continue;

          restore_bigint_from_residues(*output_residues_window, i, j, comb,
                                       residues_buffer_temp, bigint_value);
          const auto iLoc = output.LocalRow(i);
          const auto jLoc = output.LocalCol(j);
          bigint_value.to_BigFloat(output.Matrix()(iLoc, jLoc));
        }
  }

  // In a single-node case, no need to reduce anything
  if(El::mpi::Congruent(shared_memory_comm, output_comm))
    {
      // Ensure that no one will write to the window until we finish
      // (if we split Q, the window is reused multimple times)
      Scoped_Timer fence_timer(timers, "fence");
      output_residues_window->Fence();
      return;
    }

  // reduce_comm will combine processes having the same rank on a node,
  // e.g.: (0,128,256), (1,129,257),... (for 3 nodes, 128 core each)
  // TODO create communicator in constructor and reuse it?
  // Generally we don't know output_comm, in practice it's always COMM_WORLD.
  MPI_Comm_Wrapper reduce_comm_wrapper;
  auto &reduce_comm = reduce_comm_wrapper.value;
  El::mpi::Split(output_comm, shared_memory_comm.Rank(), 0, reduce_comm);
  const int reduce_rank = reduce_comm.Rank();

  // All nodes have the same number of ranks,
  // so that all reduce_comm should have the same size
  ASSERT_EQUAL(reduce_comm.Size() * shared_memory_comm.Size(),
               output_comm.Size());

  std::vector<int> global_ranks(reduce_comm.Size());
  global_ranks.at(reduce_rank) = output_comm.Rank();

  // Number of output elements owned by a given rank in reduce_comm
  std::vector<int> num_output_elements(reduce_comm.Size());
  for(int i = 0; i < height; ++i)
    for(int j = 0; j < width; ++j)
      {
        if(skip_element(i, j))
          continue;
        if(output.IsLocal(i, j))
          ++num_output_elements.at(reduce_rank);
      }
  for(int rank = 0; rank < reduce_comm.Size(); ++rank)
    {
      El::mpi::Broadcast(num_output_elements.at(rank), rank, reduce_comm);
      El::mpi::Broadcast(global_ranks.at(rank), rank, reduce_comm);
    }

  {
    Scoped_Timer reduce_timer(timers, "reduce");

    El::BigFloat bigfloat_value;
    const size_t serialized_size = bigfloat_value.SerializedSize();

    std::vector<El::byte> send_buf;
    std::vector<El::byte> recv_buf;
    // Each rank sends to (rank + offset) Q_n[i,j] for all [i,j] owned by (rank+offset)
    // TODO implement via simple MPI_Reduce?
    // Current implementation ensures uniform network load for all nodes.
    // If we switch to MPI_Reduce, we should preserve this uniformity.
    // For example, naive implementation "reduce for all ranks of node 0, then for all ranks of node 1 etc."
    // will be probably worse.
    for(int rank_offset = 1; rank_offset < reduce_comm.Size(); ++rank_offset)
      {
        Scoped_Timer iter_timer(timers,
                                "offset=" + std::to_string(rank_offset));
        const int to = (reduce_rank + rank_offset) % reduce_comm.Size();
        const int from = (reduce_comm.Size() + reduce_rank - rank_offset)
                         % reduce_comm.Size();

        // Fill send buffer
        {
          Scoped_Timer serialize_timer(timers, "serialize");
          send_buf.resize(num_output_elements.at(to) * serialized_size);
          El::byte *curr_send = send_buf.data();
          for(int i = 0; i < height; ++i)
            for(int j = 0; j < width; ++j)
              {
                if(skip_element(i, j))
                  continue;
                if(output.Owner(i, j) != global_ranks.at(to))
                  continue;

                ASSERT(curr_send + serialized_size - send_buf.data()
                         <= send_buf.size(),
                       "send buffer is too small!",
                       DEBUG_STRING(send_buf.size()),
                       DEBUG_STRING(curr_send - send_buf.data()),
                       DEBUG_STRING(serialized_size), DEBUG_STRING(i),
                       DEBUG_STRING(j), DEBUG_STRING(height),
                       DEBUG_STRING(width), DEBUG_STRING(global_ranks.at(to)));
                restore_bigint_from_residues(*output_residues_window, i, j,
                                             comb, residues_buffer_temp,
                                             bigint_value);
                bigint_value.to_BigFloat(bigfloat_value);
                bigfloat_value.Serialize(curr_send);
                curr_send += serialized_size;
              }
          ASSERT_EQUAL(curr_send - send_buf.data(), send_buf.size(),
                       "send buffer has wrong size!",
                       DEBUG_STRING(serialized_size), DEBUG_STRING(height),
                       DEBUG_STRING(width), DEBUG_STRING(global_ranks.at(to)));
        }

        // Receive buffer will receive elements for a current rank
        recv_buf.resize(num_output_elements.at(reduce_rank) * serialized_size);

        {
          Scoped_Timer mpi_sendrecv_timer(timers, "mpi_sendrecv");
          El::mpi::SendRecv(send_buf.data(), send_buf.size(), to,
                            recv_buf.data(), recv_buf.size(), from,
                            reduce_comm);
        }

        // Restore received
        {
          Scoped_Timer deserialize_timer(timers, "deserialize");
          El::byte *curr_recv = recv_buf.data();
          for(int i = 0; i < height; ++i)
            for(int j = 0; j < width; ++j)
              {
                if(skip_element(i, j))
                  continue;
                if(output.Owner(i, j) != global_ranks.at(reduce_rank))
                  continue;

                ASSERT(curr_recv - recv_buf.data() < recv_buf.size());

                bigfloat_value.Deserialize(curr_recv);

                const auto iLoc = output.LocalRow(i);
                const auto jLoc = output.LocalCol(j);
                output.Matrix()(iLoc, jLoc) += bigfloat_value;
                curr_recv += serialized_size;
              }
          ASSERT_EQUAL(curr_recv - recv_buf.data(), recv_buf.size());
        }
      }
  }

  // Ensure that no one will write to the window until we finish
  // (if we split Q, the window is reused multimple times)
  Scoped_Timer fence_timer(timers, "fence");
  output_residues_window->Fence();
}
