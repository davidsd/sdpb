#include "copy_matrix.hxx"
#include "Shared_Window_Array.hxx"
#include "assert.hxx"

// TODO: maybe it's better to store source as El::DistMatrix<El::BigFloat,STAR,STAR>?
// Check and update usages.
void copy_matrix(const El::Matrix<El::BigFloat> &source,
                 El::DistMatrix<El::BigFloat> &destination)
{
  destination.Resize(source.Height(), source.Width());
  for(int iLoc = 0; iLoc < destination.LocalHeight(); ++iLoc)
    {
      const int i = destination.GlobalRow(iLoc);
      for(int jLoc = 0; jLoc < destination.LocalWidth(); ++jLoc)
        {
          const int j = destination.GlobalCol(jLoc);
          destination.SetLocal(iLoc, jLoc, source(i, j));
        }
    }
}

void copy_matrix(const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &source,
                 El::Matrix<El::BigFloat> &destination)
{
  El::Copy(source.LockedMatrix(), destination);
}

void copy_matrix(const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &source,
                 El::DistMatrix<El::BigFloat> &destination)
{
  // NB: Calling El::Copy(source, destination) may lead to "Grids did not match" error,
  // if e.g. source is global and destination is local to some MPI group.
  copy_matrix(source.LockedMatrix(), destination);
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
      El::Copy(source, destination.Matrix());
      return;
    }
  if(comm.Rank() == 0)
    {
      for(int i = 0; i < source.Height(); ++i)
        for(int j = 0; j < source.Width(); ++j)
          destination.QueueUpdate(i, j, source.Get(i, j));
    }
  destination.ProcessQueues();
}