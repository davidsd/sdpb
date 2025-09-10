#pragma once

// Base class for Block_Matrix, Block_Vector, Block_Diagonal_Matrix etc.

#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <vector>

// Curiously Recurring Template Pattern ensures type safety for operator+= etc.
template <typename Derived> struct Abstract_Block_Matrix
{
  std::vector<El::DistMatrix<El::BigFloat>> blocks;

  virtual void add_block(size_t height, const El::Grid &grid) = 0;
  virtual ~Abstract_Block_Matrix() = default;

  // Utility functions and operators

  void set_zero()
  {
    for(auto &block : blocks)
      El::Zero(block);
  }

  // Add a constant c to each diagonal element
  void operator+=(const Derived &A)
  {
    ASSERT_EQUAL(blocks.size(), A.blocks.size());
    for(size_t b = 0; b < blocks.size(); b++)
      blocks[b] += A.blocks[b];
  }
  void operator-=(const Derived &A)
  {
    ASSERT_EQUAL(blocks.size(), A.blocks.size());
    for(size_t b = 0; b < blocks.size(); b++)
      blocks[b] -= A.blocks[b];
  }
  void operator*=(const El::BigFloat &c)
  {
    for(auto &block : blocks)
      block *= c;
  }

  // The maximal absolute value of the elements of M
  [[nodiscard]] El::BigFloat max_abs() const
  {
    El::BigFloat max = 0;
    for(auto &block : blocks)
      {
        max = std::max(El::MaxAbs(block), max);
      }
    return El::mpi::AllReduce(max, El::mpi::MAX, El::mpi::COMM_WORLD);
  }

  friend std::ostream &
  operator<<(std::ostream &os, const Abstract_Block_Matrix &A)
  {
    os << "{";
    for(auto block(A.blocks.begin()); block != A.blocks.end();)
      {
        El::Print(*block, "", ",", os);
        ++block;
        if(block != A.blocks.end())
          {
            os << ", ";
          }
      }
    os << "}";
    return os;
  }
};

// More utility functions

template <class Derived>
El::BigFloat dotu(const Abstract_Block_Matrix<Derived> &A,
                  const Abstract_Block_Matrix<Derived> &B)
{
  ASSERT_EQUAL(A.blocks.size(), B.blocks.size());
  El::BigFloat local_sum(0);

  for(size_t ii = 0; ii != A.blocks.size(); ++ii)
    {
      // FIXME: This feels slow.  It has to wait for each block
      // computation to be done before it can go to the next.
      local_sum += El::Dotu(A.blocks[ii], B.blocks[ii]);
    }
  // Make sure not to double count if blocks are distributed over more
  // than one processor.  We could also divide the sum by
  // X.blocks.front().Size().
  if(!A.blocks.empty() && A.blocks.front().Grid().Rank() != 0)
    {
      local_sum = 0;
    }
  return El::mpi::AllReduce(local_sum, El::mpi::COMM_WORLD);
}
