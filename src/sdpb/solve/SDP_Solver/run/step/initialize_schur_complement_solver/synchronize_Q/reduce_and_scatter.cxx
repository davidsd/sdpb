#include <El.hpp>

void reduce_and_scatter(const El::mpi::Comm &comm,
                        std::vector<El::byte> &send_buffer,
                        std::vector<El::byte> &receive_buffer,
                        std::vector<int> &rank_sizes)
{
  int mpi_error(
    MPI_Reduce_scatter(send_buffer.data(), receive_buffer.data(),
                       rank_sizes.data(), El::mpi::TypeMap<El::BigFloat>(),
                       El::mpi::SumOp<El::BigFloat>().op, comm.comm));

  if(mpi_error != MPI_SUCCESS)
    {
      std::vector<char> error_string(MPI_MAX_ERROR_STRING);
      int lengthOfErrorString;
      MPI_Error_string(mpi_error, error_string.data(), &lengthOfErrorString);
      El::RuntimeError(std::string(error_string.data()));
    }
}
