// Synchronize the results back to the global Q.

#include "../../../../Timers.hxx"

#include <El.hpp>

void reduce_and_scatter(const El::mpi::Comm &comm,
                        std::vector<El::byte> &sending_buffer,
                        std::vector<El::byte> &receiving_buffer,
                        std::vector<int> &rank_sizes);

void synchronize(const El::DistMatrix<El::BigFloat> &Q_group,
                 const bool &debug,
                 El::DistMatrix<El::BigFloat> &Q, Timers &timers)
{
  auto &synchronize_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.Qcomputation."
    "synchronize",debug));

  El::BigFloat zero(0);
  const size_t serialized_size(zero.SerializedSize());
  std::vector<uint8_t> serialized_zero(serialized_size);
  zero.Serialize(serialized_zero.data());

  // MPI uses 'int' for offsets.
  std::vector<int> rank_sizes(El::mpi::Size(El::mpi::COMM_WORLD));
  for(int64_t row = 0; row < Q_group.Height(); ++row)
    for(int64_t column = row; column < Q_group.Height(); ++column)
      {
        ++rank_sizes.at(Q.Owner(row, column));
      }

  std::vector<size_t> rank_offsets(rank_sizes.size() + 1, 0);
  for(size_t rank = 1; rank < rank_offsets.size(); ++rank)
    {
      rank_offsets[rank] = rank_offsets[rank - 1] + rank_sizes[rank - 1];
    }

  const size_t total_sending(rank_offsets.back() * serialized_size);
  std::vector<El::byte> sending_buffer(total_sending);

  std::vector<size_t> current_offsets(rank_offsets);
  for(int64_t row = 0; row < Q_group.Height(); ++row)
    for(int64_t column = row; column < Q_group.Height(); ++column)
      {
        size_t &offset(current_offsets.at(Q.Owner(row, column)));
        El::byte *insertion_point(sending_buffer.data()
                                  + offset * serialized_size);
        if(Q_group.IsLocal(row, column))
          {
            Q_group.GetLocal(Q_group.LocalRow(row), Q_group.LocalCol(column))
              .Serialize(insertion_point);
          }
        else
          {
            std::copy(serialized_zero.begin(), serialized_zero.end(),
                      insertion_point);
          }
        ++offset;
      }

  const size_t total_receiving(rank_sizes.at(El::mpi::Rank())
                               * serialized_size);
  std::vector<El::byte> receiving_buffer(total_receiving);

  reduce_and_scatter(El::mpi::COMM_WORLD, sending_buffer, receiving_buffer,
                     rank_sizes);

  El::byte *current_receiving(receiving_buffer.data());

  for(int64_t row = 0; row < Q.Height(); ++row)
    for(int64_t column = row; column < Q.Width(); ++column)
      {
        if(Q.Owner(row, column) == El::mpi::Rank())
          {
            El::BigFloat received;
            current_receiving = received.Deserialize(current_receiving);
            Q.SetLocal(Q.LocalRow(row), Q.LocalCol(column), received);
          }
      }
  synchronize_timer.stop();
}
