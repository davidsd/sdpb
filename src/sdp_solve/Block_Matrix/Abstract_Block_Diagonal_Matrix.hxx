#pragma once

#include "Abstract_Block_Matrix.hxx"
#include "sdpb_util/assert.hxx"

// A block-diagonal square matrix
//
//   M = Diagonal(M_0, M_1, ..., M_{bMax-1})
//
// where each block M_b is a square-matrix (of possibly different
// sizes).
template <class Derived>
struct Abstract_Block_Diagonal_Matrix : Abstract_Block_Matrix<Derived>
{
  virtual ~Abstract_Block_Diagonal_Matrix() = default;

  void add_block(const size_t height, const El::Grid &grid) override
  {
    this->blocks.emplace_back(height, height, grid);
  }

  // Add a constant c to each diagonal element
  void add_diagonal(const El::BigFloat &c)
  {
    for(auto &block : this->blocks)
      ShiftDiagonal(block, c);
  }
  void symmetrize()
  {
    for(auto &block : this->blocks)
      {
        // FIXME: This feels expensive

        // We can not use El::MakeSymmetric() because that just copies
        // the lower part to the upper part.  We need to average the
        // upper and lower parts.
        block *= 0.5;
        El::DistMatrix<El::BigFloat> transpose(block.Grid());
        El::Transpose(block, transpose, false);
        block += transpose;
      }
  }

  [[nodiscard]] El::BigFloat trace() const
  {
    El::BigFloat result = 0;
    for(auto &block : this->blocks)
      {
        const auto block_trace = El::Trace(block);
        // to avoid double-counting, we accumulate results on rank 0 of each block
        if(block.DistComm().Rank() == 0)
          result += block_trace;
      }
    return El::mpi::AllReduce(result, El::mpi::COMM_WORLD);
  }
};

// C := alpha*A*B + beta*C
template <class Derived>
void scale_multiply_add(const El::BigFloat &alpha,
                        const Abstract_Block_Diagonal_Matrix<Derived> &A,
                        const Abstract_Block_Diagonal_Matrix<Derived> &B,
                        const El::BigFloat &beta,
                        Abstract_Block_Diagonal_Matrix<Derived> &C)
{
  ASSERT_EQUAL(A.blocks.size(), B.blocks.size());
  ASSERT_EQUAL(A.blocks.size(), C.blocks.size());
  for(size_t block = 0; block < A.blocks.size(); ++block)
    {
      El::Gemm(El::OrientationNS::NORMAL, El::OrientationNS::NORMAL, alpha,
               A.blocks[block], B.blocks[block], beta, C.blocks[block]);
    }
}

// C := A B
template <class Derived>
void multiply(const Abstract_Block_Diagonal_Matrix<Derived> &A,
              const Abstract_Block_Diagonal_Matrix<Derived> &B,
              Abstract_Block_Diagonal_Matrix<Derived> &C)
{
  scale_multiply_add(El::BigFloat(1), A, B, El::BigFloat(0), C);
}

// X := ACholesky^{-T} ACholesky^{-1} X = A^{-1} X
template <class Derived>
void cholesky_solve(const Abstract_Block_Diagonal_Matrix<Derived> &ACholesky,
                    Abstract_Block_Diagonal_Matrix<Derived> &X)
{
  ASSERT_EQUAL(ACholesky.blocks.size(), X.blocks.size());
  for(size_t b = 0; b < X.blocks.size(); b++)
    {
      El::cholesky::SolveAfter(El::UpperOrLowerNS::LOWER,
                               El::OrientationNS::NORMAL, ACholesky.blocks[b],
                               X.blocks[b]);
    }
}

// B := L^{-1} B,
// where L is the result of a previous cholesky factorization.
// Note that this is different from computing the solution to A B == (L L^T) B
template <class Derived_L, class Derived_B>
void lower_triangular_solve(
  const Abstract_Block_Diagonal_Matrix<Derived_L> &L_cholesky,
  Abstract_Block_Matrix<Derived_B> &B)
{
  ASSERT_EQUAL(L_cholesky.blocks.size(), B.blocks.size());
  for(size_t block = 0; block < L_cholesky.blocks.size(); block++)
    {
      El::Trsm(El::LeftOrRightNS::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::NORMAL, El::UnitOrNonUnitNS::NON_UNIT,
               El::BigFloat(1), L_cholesky.blocks[block], B.blocks[block]);
    }
}

// v := L^{-T} v, where L is lower-triangular
template <class Derived_L, class Derived_B>
void lower_triangular_transpose_solve(const Abstract_Block_Diagonal_Matrix<Derived_L> &L,
                                      Abstract_Block_Matrix<Derived_B> &v)
{
  for(size_t b = 0; b < L.blocks.size(); b++)
    {
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::OrientationNS::TRANSPOSE, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L.blocks[b], v.blocks[b]);
    }
}

// A := L^{-1} A L^{-T}
template <class Derived>
void lower_triangular_inverse_congruence(
  const Abstract_Block_Diagonal_Matrix<Derived> &L,
  Abstract_Block_Diagonal_Matrix<Derived> &A)
{
  for(size_t b = 0; b < A.blocks.size(); b++)
    {
      El::Trsm(El::LeftOrRight::RIGHT, El::UpperOrLowerNS::LOWER,
               El::Orientation::TRANSPOSE, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L.blocks[b], A.blocks[b]);
      El::Trsm(El::LeftOrRight::LEFT, El::UpperOrLowerNS::LOWER,
               El::Orientation::NORMAL, El::UnitOrNonUnit::NON_UNIT,
               El::BigFloat(1), L.blocks[b], A.blocks[b]);
    }
}

// Tr(A B), where A and B are symmetric
template <class Derived>
El::BigFloat
frobenius_product_symmetric(const Abstract_Block_Diagonal_Matrix<Derived> &A,
                            const Abstract_Block_Diagonal_Matrix<Derived> &B)
{
  return dotu(A, B);
}

// (X + dX) . (Y + dY), where X, dX, Y, dY are symmetric
// BlockDiagonalMatrices and '.' is the Frobenius product.
template <class Derived>
El::BigFloat
frobenius_product_of_sums(const Abstract_Block_Diagonal_Matrix<Derived> &X,
                          const Abstract_Block_Diagonal_Matrix<Derived> &dX,
                          const Abstract_Block_Diagonal_Matrix<Derived> &Y,
                          const Abstract_Block_Diagonal_Matrix<Derived> &dY)
{
  ASSERT_EQUAL(X.blocks.size(), dX.blocks.size());
  ASSERT_EQUAL(X.blocks.size(), Y.blocks.size());
  ASSERT_EQUAL(Y.blocks.size(), dY.blocks.size());
  El::BigFloat local_sum(0);
  for(size_t b = 0; b < X.blocks.size(); b++)
    {
      // FIXME: This can be sped up by not have intermediate results.
      // It may require looking into the implementation of Dotu.
      El::DistMatrix<El::BigFloat> X_dX(X.blocks[b]);
      X_dX += dX.blocks[b];
      El::DistMatrix<El::BigFloat> Y_dY(Y.blocks[b]);
      Y_dY += dY.blocks[b];
      local_sum += Dotu(X_dX, Y_dY);
    }

  // Make sure not to double count if blocks are distributed over more
  // than one processor.  We could also divide the sum by
  // X.blocks.front().Size().
  if(!X.blocks.empty() && X.blocks.front().Grid().Rank() != 0)
    {
      local_sum = 0;
    }
  return El::mpi::AllReduce(local_sum, El::mpi::COMM_WORLD);
}
