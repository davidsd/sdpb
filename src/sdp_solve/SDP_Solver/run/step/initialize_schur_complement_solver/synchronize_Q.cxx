// Synchronize the results back to the global Q.

#include "../../../../../Timers.hxx"

#include <El.hpp>

#include <algorithm>

namespace
{
  void check_mpi_error(const int &mpi_error)
  {
    if(mpi_error != MPI_SUCCESS)
      {
        std::vector<char> error_string(MPI_MAX_ERROR_STRING);
        int lengthOfErrorString;
        MPI_Error_string(mpi_error, error_string.data(), &lengthOfErrorString);
        El::RuntimeError(std::string(error_string.data()));
      }
  }
}

void synchronize_Q(El::DistMatrix<El::BigFloat> &Q,
                   const El::DistMatrix<El::BigFloat> &Q_group, Timers &timers)
{
  Scoped_Timer timer(timers, "synchronize_Q");

  const int total_ranks(El::mpi::Size(El::mpi::COMM_WORLD));
  // Special case serial case
  if(total_ranks == 1)
    {
      for(int64_t row = 0; row < Q_group.Height(); ++row)
        for(int64_t column = row; column < Q_group.Height(); ++column)
          {
            Q.SetLocal(Q.LocalRow(row), Q.LocalCol(column),
                       Q_group.GetLocal(Q_group.LocalRow(row),
                                        Q_group.LocalCol(column)));
          }
      return;
    }

  El::BigFloat zero(0);
  size_t serialized_size(zero.SerializedSize());
  std::vector<uint8_t> serialized_zero(serialized_size);
  zero.Serialize(serialized_zero.data());

  // This is an re-implementation of MPI_Reduce_scatter
  // using the ring algorithm as found in OpenMPI.

  // We re-implement MPI_Reduce_scatter because we can get away with
  // significantly less memory use by not constructing the full send
  // buffer beforehand.  Also, for large blocks, we can skip some
  // elements when summing because those processors do not have
  // contributions for all of Q.

  // MPI uses 'int' for message sizes.
  std::vector<int> rank_sizes(total_ranks);
  for(int64_t row = 0; row < Q_group.Height(); ++row)
    for(int64_t column = row; column < Q_group.Height(); ++column)
      {
        ++rank_sizes.at(Q.Owner(row, column));
      }

  int max_buffer_size(*std::max_element(rank_sizes.begin(), rank_sizes.end())
                      * serialized_size);
  std::array<std::vector<El::byte>, 2> send_receive(
    {std::vector<El::byte>(max_buffer_size),
     std::vector<El::byte>(max_buffer_size)});

  const int rank(El::mpi::Rank(El::mpi::COMM_WORLD));
  const int send_to_rank((rank + 1) % total_ranks),
    receive_from_rank((total_ranks + rank - 1) % total_ranks);

  // Initial async receive
  int final_receive_destination((total_ranks + rank - 2) % total_ranks);
  std::array<MPI_Request, 2> receive_requests;
  check_mpi_error(MPI_Irecv(send_receive[0].data(),
                            rank_sizes[final_receive_destination],
                            El::mpi::TypeMap<El::BigFloat>(),
                            receive_from_rank, final_receive_destination,
                            El::mpi::COMM_WORLD.comm, &receive_requests[0]));

  // Initial fill of send buffer
  int final_send_destination((total_ranks + rank - 1) % total_ranks);

  {
    El::byte *insertion_point(send_receive[1].data());
    for(int64_t row = 0; row < Q_group.Height(); ++row)
      for(int64_t column = row; column < Q_group.Height(); ++column)
        {
          if(Q.Owner(row, column) == final_send_destination)
            {
              if(Q_group.IsLocal(row, column))
                {
                  Q_group
                    .GetLocal(Q_group.LocalRow(row), Q_group.LocalCol(column))
                    .Serialize(insertion_point);
                }
              else
                {
                  std::copy(serialized_zero.begin(), serialized_zero.end(),
                            insertion_point);
                }
              insertion_point += serialized_size;
            }
        }
  }

  check_mpi_error(MPI_Send(send_receive[1].data(),
                           rank_sizes[final_send_destination],
                           El::mpi::TypeMap<El::BigFloat>(), send_to_rank,
                           final_send_destination, El::mpi::COMM_WORLD.comm));

  // Loop over all remaining intermediate ranks
  for(int rank_offset(2); rank_offset < total_ranks; ++rank_offset)
    {
      {
        auto &receive_buffer(send_receive[(rank_offset + 1) % 2]);

        final_receive_destination
          = (total_ranks + rank - (rank_offset + 1)) % total_ranks;

        check_mpi_error(MPI_Irecv(
          receive_buffer.data(), rank_sizes[final_receive_destination],
          El::mpi::TypeMap<El::BigFloat>(), receive_from_rank,
          final_receive_destination, El::mpi::COMM_WORLD.comm,
          &receive_requests[(rank_offset + 1) % 2]));
      }
      {
        // This waits for the receive from a previous iteration, not the
        // one we just initiated.

        // We do not cancel sends, so no need to check status.
        check_mpi_error(
          MPI_Wait(&receive_requests[rank_offset % 2], MPI_STATUS_IGNORE));

        final_send_destination
          = (total_ranks + rank - rank_offset) % total_ranks;
        auto &receive_then_send_buffer(send_receive[rank_offset % 2]);
        El::byte *current_receiving(receive_then_send_buffer.data());
        for(int64_t row = 0; row < Q_group.Height(); ++row)
          for(int64_t column = row; column < Q_group.Height(); ++column)
            {
              if(Q.Owner(row, column) == final_send_destination)
                {
                  if(Q_group.IsLocal(row, column))
                    {
                      El::BigFloat received;
                      received.Deserialize(current_receiving);
                      received += Q_group.GetLocal(Q_group.LocalRow(row),
                                                   Q_group.LocalCol(column));
                      received.Serialize(current_receiving);
                    }
                  current_receiving += serialized_size;
                }
            }
        check_mpi_error(MPI_Send(
          receive_then_send_buffer.data(), rank_sizes[final_send_destination],
          El::mpi::TypeMap<El::BigFloat>(), send_to_rank,
          final_send_destination, El::mpi::COMM_WORLD.comm));
      }
    }
  // Add the local contribution to the last message and put it into
  // the global Q.

  check_mpi_error(
    MPI_Wait(&receive_requests[total_ranks % 2], MPI_STATUS_IGNORE));
  El::byte *current_receiving(send_receive[total_ranks % 2].data());
  for(int64_t row = 0; row < Q_group.Height(); ++row)
    for(int64_t column = row; column < Q_group.Height(); ++column)
      {
        if(Q.Owner(row, column) == rank)
          {
            El::BigFloat received;
            received.Deserialize(current_receiving);
            if(Q_group.IsLocal(row, column))
              {
                received += Q_group.GetLocal(Q_group.LocalRow(row),
                                             Q_group.LocalCol(column));
              }
            Q.SetLocal(Q.LocalRow(row), Q.LocalCol(column), received);
            current_receiving += serialized_size;
          }
      }
}
