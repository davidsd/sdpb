#pragma once

#include <El.hpp>

// source: same Matrix copied over all ranks of destination.DistComm()
// destination: DistMatrix
void copy_matrix(const El::Matrix<El::BigFloat> &source,
                 El::DistMatrix<El::BigFloat> &destination);

void copy_matrix(const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &source,
                 El::Matrix<El::BigFloat> &destination);

void copy_matrix(const El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &source,
                 El::DistMatrix<El::BigFloat> &destination);

// source: Matrix initialized at comm.Rank() = 0
// destination: DistMatrix of the same size, elements will be copied from comm root.
// comm should be equal to destination.DistComm()
// (this argument is added for safety check)
void copy_matrix_from_root(const El::Matrix<El::BigFloat> &source,
                           El::DistMatrix<El::BigFloat> &destination,
                           const El::mpi::Comm &comm);
