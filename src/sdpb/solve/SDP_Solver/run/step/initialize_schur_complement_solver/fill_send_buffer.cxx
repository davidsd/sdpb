// Synchronize the results back to the global Q.

#include "../../../../../../Timers.hxx"

#include <El.hpp>

void fill_send_buffer(const El::DistMatrix<El::BigFloat> &Q,
                      const El::DistMatrix<El::BigFloat> &Q_group,
                      std::vector<El::byte> &send_buffer,
                      std::vector<int> &rank_sizes, size_t &serialized_size,
                      Timers &timers)
{
  auto &fill_send_buffers_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.Q.fill_send_buffers"));

  El::BigFloat zero(0);
  serialized_size = zero.SerializedSize();
  std::vector<uint8_t> serialized_zero(serialized_size);
  zero.Serialize(serialized_zero.data());

  // MPI uses 'int' for offsets.
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
  send_buffer.resize(total_sending);

  std::vector<size_t> current_offsets(rank_offsets);
  for(int64_t row = 0; row < Q_group.Height(); ++row)
    for(int64_t column = row; column < Q_group.Height(); ++column)
      {
        size_t &offset(current_offsets.at(Q.Owner(row, column)));
        El::byte *insertion_point(send_buffer.data()
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
  fill_send_buffers_timer.stop();
}
