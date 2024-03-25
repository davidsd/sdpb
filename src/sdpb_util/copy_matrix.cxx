#include "copy_matrix.hxx"
#include "Shared_Window_Array.hxx"
#include "assert.hxx"

void copy_matrix(const El::Matrix<El::BigFloat> &source,
                 El::DistMatrix<El::BigFloat> &destination)
{
  destination.Resize(source.Height(), source.Width());
  for(int64_t row(0); row < destination.LocalHeight(); ++row)
    {
      int64_t global_row(destination.GlobalRow(row));
      for(int64_t column(0); column < destination.LocalWidth(); ++column)
        {
          int64_t global_column(destination.GlobalCol(column));
          destination.SetLocal(row, column, source(global_row, global_column));
        }
    }
}

void copy_matrix(const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &source,
                 El::Matrix<El::BigFloat> &destination)
{
  destination.Resize(source.LocalHeight(), source.LocalWidth());
  for(int64_t row(0); row < source.LocalHeight(); ++row)
    {
      for(int64_t column(0); column < source.LocalWidth(); ++column)
        {
          destination(row, column) = source.GetLocal(row, column);
        }
    }
}

void copy_matrix(const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &source,
                 El::DistMatrix<El::BigFloat> &destination)
{
  destination.Resize(source.LocalHeight(), source.LocalWidth());
  for(int64_t row(0); row < destination.LocalHeight(); ++row)
    {
      const int64_t global_row(destination.GlobalRow(row));
      for(int64_t column(0); column < destination.LocalWidth(); ++column)
        {
          const int64_t global_column(destination.GlobalCol(column));
          destination.SetLocal(row, column,
                               source.GetLocal(global_row, global_column));
        }
    }
}

// source: Matrix initialized at comm.Rank() = 0
// destination: DistMatrix of the same size, elements will be copied from comm root.
// comm should be equal to output.DistComm()
// (this argument is added for safety check)
//
// In our real-world tests, this version is ~3x slower than copy_matrix_from_root_impl_shared_window()
void copy_matrix_from_root_impl_send_recv(
  const El::Matrix<El::BigFloat> &source,
  El::DistMatrix<El::BigFloat> &destination, const El::mpi::Comm &comm)
{
  const El::Int root = 0;
  const El::Int rank = comm.Rank();

  for(El::Int global_row = 0; global_row < destination.Height(); ++global_row)
    for(El::Int global_col = 0; global_col < destination.Width(); ++global_col)
      {
        const auto from = root;
        const auto to = destination.Owner(global_row, global_col);
        // special case: copy from root to root, no MPI messaging needed
        if(to == from)
          {
            if(rank == to)
              {
                destination.Set(global_row, global_col,
                                source.Get(global_row, global_col));
              }
            continue;
          }
        // send from root
        if(rank == from)
          {
            El::mpi::Send(source.Get(global_row, global_col), to, comm);
          }
        // recieve by element owner
        if(rank == to)
          {
            destination.Set(global_row, global_col,
                            El::mpi::Recv<El::BigFloat>(from, comm));
          }
      }
}

// source: Matrix initialized at comm.Rank() = 0
// destination: DistMatrix of the same size, elements will be copied from comm root.
// comm should be equal to output.DistComm()
// (this argument is added for safety check)
void copy_matrix_from_root_impl_shared_window(
  const El::Matrix<El::BigFloat> &source,
  El::DistMatrix<El::BigFloat> &destination, const El::mpi::Comm &comm)
{
  const El::Int root = 0;
  const El::Int rank = comm.Rank();

  const auto bigfloat_size_bytes
    = El::BigFloat(1).SerializedSize() / sizeof(El::byte);

  Shared_Window_Array<El::byte> window(
    comm, destination.Height() * destination.Width() * bigfloat_size_bytes);

  const int ldim = destination.Height();
  auto get_buffer = [&](int i, int j) -> El::byte * {
    return window.data + (i + j * ldim) * bigfloat_size_bytes;
  };

  // Serialize input matrix to shared memory window
  if(rank == root)
    {
      for(int i = 0; i < source.Height(); ++i)
        for(int j = 0; j < source.Width(); ++j)
          {
            source.Get(i, j).Serialize(get_buffer(i, j));
          }
    }
  window.Fence();

  // Deserialize to output
  for(int iLoc = 0; iLoc < destination.LocalHeight(); ++iLoc)
    for(int jLoc = 0; jLoc < destination.LocalWidth(); ++jLoc)
      {
        int i = destination.GlobalRow(iLoc);
        int j = destination.GlobalCol(jLoc);
        if(rank == root)
          {
            destination.SetLocal(iLoc, jLoc, source.Get(i, j));
          }
        else
          {
            El::BigFloat value;
            value.Deserialize(get_buffer(i, j));
            destination.SetLocal(iLoc, jLoc, value);
          }
      }
}

// source: Matrix initialized at comm.Rank() = 0
// destination: DistMatrix of the same size, elements will be copied from comm root.
// comm should be equal to output.DistComm()
// (this argument is added for safety check)
void copy_matrix_from_root(const El::Matrix<El::BigFloat> &source,
                           El::DistMatrix<El::BigFloat> &destination,
                           const El::mpi::Comm &comm)
{
  ASSERT(El::mpi::Congruent(comm, destination.DistComm()),
         "Wrong communicator for copy_matrix_from_root(); "
         "use output.DistComm()");

  if(comm.Rank() == 0)
    {
      // TODO set matrix size instead?
      ASSERT(source.Height() == destination.Height()
               && source.Width() == destination.Width(),
             "Incompatible matrix sizes: ", El::DimsString(source, "source"),
             ", ", El::DimsString(destination, "destination"));
    }

  if(comm.Size() == 1)
    {
      copy_matrix(source, destination);
      return;
    }

  // In our real-world tests, copy_matrix_from_root_impl_shared_window() was ~3x faster than copy_matrix_from_root_impl_send_recv()
  // Thus we use it when possible, i.e. when all ranks are on the same node

  El::mpi::Comm shared_memory_comm;
  MPI_Comm_split_type(comm.comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                      &shared_memory_comm.comm);

  if(El::mpi::Congruent(comm, shared_memory_comm))
    copy_matrix_from_root_impl_shared_window(source, destination, comm);
  else
    copy_matrix_from_root_impl_send_recv(source, destination, comm);
  El::mpi::Free(shared_memory_comm);
}