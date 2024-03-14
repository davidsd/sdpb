#pragma once

#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <El.hpp>

#include <algorithm>
#include <optional>

inline void check_mpi_error(const int &mpi_error)
{
  if(mpi_error != MPI_SUCCESS)
    {
      std::vector<char> error_string(MPI_MAX_ERROR_STRING);
      int lengthOfErrorString;
      MPI_Error_string(mpi_error, error_string.data(), &lengthOfErrorString);
      RUNTIME_ERROR(std::string(error_string.data()));
    }
}

// input is a matrix distributed over all ranks on a node.
// output is a global matrix.
// We set output = (sum of input from all nodes).
inline void
reduce_scatter(El::DistMatrix<El::BigFloat> &output,
               const El::DistMatrix<El::BigFloat> &input, Timers &timers,
               const std::optional<El::UpperOrLower> &uplo = std::nullopt)
{
  Scoped_Timer timer(timers, "reduce_scatter");

  ASSERT_EQUAL(input.Height(), output.Height());
  ASSERT_EQUAL(input.Width(), output.Width());

  if(input.DistComm().Size() == output.DistComm().Size())
    ASSERT(El::mpi::Congruent(input.DistComm(), output.DistComm()));
  else
    ASSERT(input.DistComm().Size() < output.DistComm().Size());

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
  const int total_ranks = output_comm.Size();
  // Special case: input and output are spread over the same ranks,
  // no need for reduce-scatter
  if(input.DistComm().Size() == output_comm.Size())
    {
      for(El::Int row = 0; row < input.Height(); ++row)
        for(El::Int column = 0; column < input.Width(); ++column)
          {
            if(skip_element(row, column))
              continue;
            output.SetLocal(
              output.LocalRow(row), output.LocalCol(column),
              input.GetLocal(input.LocalRow(row), input.LocalCol(column)));
          }
      return;
    }

  const El::BigFloat zero(0);
  const size_t serialized_size = zero.SerializedSize();
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
  for(El::Int row = 0; row < input.Height(); ++row)
    for(El::Int column = 0; column < input.Width(); ++column)
      {
        if(skip_element(row, column))
          continue;
        ++rank_sizes.at(output.Owner(row, column));
      }

  const int max_buffer_size
    = *std::max_element(rank_sizes.begin(), rank_sizes.end())
      * serialized_size;
  std::array<std::vector<El::byte>, 2> send_receive(
    {std::vector<El::byte>(max_buffer_size),
     std::vector<El::byte>(max_buffer_size)});

  const int rank = output_comm.Rank();
  const int send_to_rank = (rank + 1) % total_ranks;
  const int receive_from_rank = (total_ranks + rank - 1) % total_ranks;

  // Initial async receive
  int final_receive_destination = (total_ranks + rank - 2) % total_ranks;
  std::array<MPI_Request, 2> receive_requests{};
  check_mpi_error(MPI_Irecv(
    send_receive[0].data(), rank_sizes[final_receive_destination],
    El::mpi::TypeMap<El::BigFloat>(), receive_from_rank,
    final_receive_destination, output_comm.comm, &receive_requests[0]));

  // Initial fill of send buffer
  int final_send_destination = (total_ranks + rank - 1) % total_ranks;

  {
    El::byte *insertion_point(send_receive[1].data());
    for(El::Int row = 0; row < input.Height(); ++row)
      for(El::Int column = 0; column < input.Width(); ++column)
        {
          if(skip_element(row, column))
            continue;
          if(output.Owner(row, column) == final_send_destination)
            {
              if(input.IsLocal(row, column))
                {
                  input.GetLocal(input.LocalRow(row), input.LocalCol(column))
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
                           final_send_destination, output_comm.comm));

  // Loop over all remaining intermediate ranks
  for(int rank_offset(2); rank_offset < total_ranks; ++rank_offset)
    {
      {
        auto &receive_buffer = send_receive[(rank_offset + 1) % 2];

        final_receive_destination
          = (total_ranks + rank - (rank_offset + 1)) % total_ranks;

        check_mpi_error(MPI_Irecv(
          receive_buffer.data(), rank_sizes[final_receive_destination],
          El::mpi::TypeMap<El::BigFloat>(), receive_from_rank,
          final_receive_destination, output_comm.comm,
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
        auto &receive_then_send_buffer = send_receive[rank_offset % 2];
        El::byte *current_receiving = receive_then_send_buffer.data();
        for(El::Int row = 0; row < input.Height(); ++row)
          for(El::Int column = 0; column < input.Width(); ++column)
            {
              if(skip_element(row, column))
                continue;
              if(output.Owner(row, column) == final_send_destination)
                {
                  if(input.IsLocal(row, column))
                    {
                      El::BigFloat received;
                      received.Deserialize(current_receiving);
                      received += input.GetLocal(input.LocalRow(row),
                                                 input.LocalCol(column));
                      received.Serialize(current_receiving);
                    }
                  current_receiving += serialized_size;
                }
            }
        check_mpi_error(MPI_Send(
          receive_then_send_buffer.data(), rank_sizes[final_send_destination],
          El::mpi::TypeMap<El::BigFloat>(), send_to_rank,
          final_send_destination, output_comm.comm));
      }
    }
  // Add the local contribution to the last message and put it into
  // the global Q.

  check_mpi_error(
    MPI_Wait(&receive_requests[total_ranks % 2], MPI_STATUS_IGNORE));
  El::byte *current_receiving = send_receive[total_ranks % 2].data();
  for(El::Int row = 0; row < input.Height(); ++row)
    for(El::Int column = 0; column < input.Width(); ++column)
      {
        if(skip_element(row, column))
          continue;
        if(output.Owner(row, column) == rank)
          {
            El::BigFloat received;
            received.Deserialize(current_receiving);
            if(input.IsLocal(row, column))
              {
                received += input.GetLocal(input.LocalRow(row),
                                           input.LocalCol(column));
              }
            output.SetLocal(output.LocalRow(row), output.LocalCol(column),
                            received);
            current_receiving += serialized_size;
          }
      }
}
