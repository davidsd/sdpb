// Synchronize the results back to the global Q.

#include "../../../../../../../Timers.hxx"

#include <El.hpp>

void reduce_and_scatter(const El::mpi::Comm &comm,
                        std::vector<El::byte> &send_buffer,
                        std::vector<El::byte> &receive_buffer,
                        std::vector<int> &rank_sizes);

void synchronize_Q(std::vector<El::byte> &send_buffer,
                   std::vector<int> &rank_sizes,
                   const size_t &serialized_size,
                   El::DistMatrix<El::BigFloat> &Q, Timers &timers)
{
  auto &synchronize_Q_timer(timers.add_and_start(
    "run.step.initializeSchurComplementSolver.Qcomputation."
    "synchronize_Q"));
  const int rank(El::mpi::Rank());
  const size_t total_receiving(rank_sizes.at(rank) * serialized_size);
  std::vector<El::byte> receive_buffer(total_receiving);

  reduce_and_scatter(El::mpi::COMM_WORLD, send_buffer, receive_buffer,
                     rank_sizes);

  El::byte *current_receiving(receive_buffer.data());

  for(int64_t row = 0; row < Q.Height(); ++row)
    for(int64_t column = row; column < Q.Width(); ++column)
      {
        if(Q.Owner(row, column) == rank)
          {
            El::BigFloat received;
            current_receiving = received.Deserialize(current_receiving);
            Q.SetLocal(Q.LocalRow(row), Q.LocalCol(column), received);
          }
      }
  synchronize_Q_timer.stop();
}
