#pragma once

#include "sdp_solve/Block_Matrix/Block_Diagonal_Matrix.hxx"
#include "sdp_solve/Block_Matrix/Block_Matrix.hxx"

#include <type_traits>

// Helper classes and functions to convert BigFloat matrices to big integer matrices
// by normalizing each row or column and multiplyng each element by 2^bits (bits = El::gmp::Precision()).
// This is used for bigint arithmetics (fmpz_t, see FLINT library)
// e.g. for calculating Q := P^T * P via bigint arithmetics Chinese remainder theorem and BLAS.
// Usage:
// Block_Matrix P; DistMatrix Q;
// const Block_Matrix_Normalizer normalizer(P, El::gmp::Precision(), El::mpi::COMM_WORLD);
// normalizer.normalize_and_shift(P);
// ... // calculate Q := P^T * P
// restore_syrk_result(uplo, normalizer, Q);

// Helper for static_assert
// TODO move to sdpb_util
template <typename T = void> inline constexpr bool always_false_v = false;

// Matrix_Normalization_Kind

enum class Matrix_Normalization_Kind
{
  // Normalize each matrix column
  COLUMNS,
  // Normalize each matrix row
  ROWS
};
inline std::ostream &
operator<<(std::ostream &os, const Matrix_Normalization_Kind kind)
{
  switch(kind)
    {
    case Matrix_Normalization_Kind::COLUMNS: return os << "COLUMNS";
    case Matrix_Normalization_Kind::ROWS: return os << "ROWS";
    default: LOGIC_ERROR("Unknown Matrix_Normalization_Kind");
    }
}

inline Matrix_Normalization_Kind
transpose(const El::Orientation orientation,
          const Matrix_Normalization_Kind kind)
{
  switch(orientation)
    {
    case El::NORMAL: return kind;
    case El::TRANSPOSE:
      case El::ADJOINT: {
        const int as_int = static_cast<int>(kind);
        ASSERT(as_int == 0 || as_int == 1);
        return static_cast<Matrix_Normalization_Kind>(1 - as_int);
      }
    default: LOGIC_ERROR("Unknown orientation");
    }
}

// Helper functions

namespace Matrix_Normalizer_Util
{
  template <class T> int local_height(const El::Matrix<T> &m)
  {
    return m.Height();
  }
  template <class T> int local_width(const El::Matrix<T> &m)
  {
    return m.Width();
  }
  template <class T> int local_height(const El::AbstractDistMatrix<T> &m)
  {
    return m.LocalHeight();
  }
  template <class T> int local_width(const El::AbstractDistMatrix<T> &m)
  {
    return m.LocalWidth();
  }
  template <class TMatrix> int global_height(const TMatrix &m)
  {
    return m.Height();
  }
  template <class TMatrix> int global_width(const TMatrix &m)
  {
    return m.Width();
  }

  template <class T> int global_row(const El::Matrix<T> &, const int iLoc)
  {
    return iLoc;
  }
  template <class T>
  int global_row(const El::AbstractDistMatrix<T> &matrix, int iLoc)
  {
    return matrix.GlobalRow(iLoc);
  }

  template <class T> int global_col(const El::Matrix<T> &, const int jLoc)
  {
    return jLoc;
  }
  template <class T>
  int global_col(const El::AbstractDistMatrix<T> &matrix, int jLoc)
  {
    return matrix.GlobalCol(jLoc);
  }
  inline int global_col(const Block_Matrix &matrix, const int jLoc)
  {
    ASSERT(!matrix.blocks.empty());
    return global_col(matrix.blocks.at(0), jLoc);
  }
  template <class T>
  T &get_local_ref(El::Matrix<T> &matrix, int iLoc, int jLoc)
  {
    return matrix.Ref(iLoc, jLoc);
  }
  template <class T>
  T &get_local_ref(El::AbstractDistMatrix<T> &matrix, int iLoc, int jLoc)
  {
    return get_local_ref(matrix.Matrix(), iLoc, jLoc);
  }

  template <class T>
  const T &get_local_cref(const El::Matrix<T> &matrix, int iLoc, int jLoc)
  {
    return matrix.CRef(iLoc, jLoc);
  }
  template <class T>
  const T &
  get_local_cref(const El::AbstractDistMatrix<T> &matrix, int iLoc, int jLoc)
  {
    return matrix.GetLocalCRef(iLoc, jLoc);
  }

  template <class T>
  void set_local(El::Matrix<T> &matrix, int iLoc, int jLoc, const T &value)
  {
    matrix.Set(iLoc, jLoc, value);
  }
  template <class T>
  void set_local(El::AbstractDistMatrix<T> &matrix, int iLoc, int jLoc,
                 const T &value)
  {
    matrix.SetLocal(iLoc, jLoc, value);
  }

  template <
    class TMatrix,
    typename = std::enable_if_t<
      std::is_same_v<El::Matrix<El::BigFloat>, TMatrix>
      || std::is_base_of_v<El::AbstractDistMatrix<El::BigFloat>, TMatrix>>>
  void update_local_elements(
    TMatrix &matrix,
    const std::function<void(int, int, El::BigFloat &)> update_local_element,
    const std::optional<El::UpperOrLower> &uplo = std::nullopt)
  {
    for(int iLoc = 0; iLoc < local_height(matrix); ++iLoc)
      {
        for(int jLoc = 0; jLoc < local_width(matrix); ++jLoc)
          {
            if(uplo.has_value())
              {
                const int i = global_row(matrix, iLoc);
                const int j = global_col(matrix, jLoc);
                if(uplo == El::UPPER && i > j)
                  continue;
                if(uplo == El::LOWER && i < j)
                  continue;
              }
            update_local_element(iLoc, jLoc,
                                 get_local_ref(matrix, iLoc, jLoc));
          }
      }
  }

  template <class TDerived>
  void update_local_elements(
    Abstract_Block_Matrix<TDerived> &block_matrix,
    const std::function<void(int, int, El::BigFloat &)> update_local_element,
    const std::optional<El::UpperOrLower> &uplo = std::nullopt)
  {
    for(auto &matrix : block_matrix.blocks)
      update_local_elements(matrix, update_local_element, uplo);
  }

  // Size of norms vector
  template <Matrix_Normalization_Kind kind, class TMatrix>
  int norms_size(const TMatrix &matrix)
  {
    if constexpr(kind == Matrix_Normalization_Kind::COLUMNS)
      return global_width(matrix);
    else if constexpr(kind == Matrix_Normalization_Kind::ROWS)
      return global_height(matrix);
    else
      static_assert(always_false_v<>, "Unknown Matrix_Normalization_Kind");
    return -1;
  }
  template <Matrix_Normalization_Kind kind, class TMatrix>
  int global_row_or_col(const TMatrix &matrix, const int iLoc, const int jLoc)
  {
    if constexpr(kind == Matrix_Normalization_Kind::COLUMNS)
      return global_col(matrix, jLoc);
    else if constexpr(kind == Matrix_Normalization_Kind::ROWS)
      return global_row(matrix, iLoc);
    else
      static_assert(always_false_v<>, "Unknown Matrix_Normalization_Kind");
    return -1;
  }

  // For each column:
  // Sum of squares of elements stored in this rank
  template <Matrix_Normalization_Kind kind, class TMatrix>
  void add_local_norms_squared(const TMatrix &matrix,
                               const std::optional<El::UpperOrLower> &uplo,
                               std::vector<El::BigFloat> &local_norms_squared)
  {
    ASSERT_EQUAL(local_norms_squared.size(), norms_size<kind>(matrix),
                 DEBUG_STRING(matrix.Height()), DEBUG_STRING(matrix.Width()));
    for(int iLoc = 0; iLoc < local_height(matrix); ++iLoc)
      for(int jLoc = 0; jLoc < local_width(matrix); ++jLoc)
        {
          if(uplo.has_value())
            {
              const int i = global_row(matrix, iLoc);
              const int j = global_col(matrix, jLoc);
              if(uplo == El::UPPER && i > j)
                continue;
              if(uplo == El::LOWER && i < j)
                continue;
            }
          const auto &value = get_local_cref(matrix, iLoc, jLoc);
          auto row_or_col = global_row_or_col<kind>(matrix, iLoc, jLoc);
          local_norms_squared.at(row_or_col) += value * value;
        }
  }

  inline void sqrt(std::vector<El::BigFloat> &values)
  {
    for(auto &value : values)
      value = El::Sqrt(value);
  }

  template <Matrix_Normalization_Kind kind, class TMatrix>
  std::vector<El::BigFloat>
  compute_norms(const TMatrix &matrix,
                const std::optional<El::UpperOrLower> &uplo = std::nullopt)
  {
    if constexpr(std::is_same_v<El::Matrix<El::BigFloat>, TMatrix>)
      {
        std::vector<El::BigFloat> norms(norms_size<kind>(matrix), 0);
        add_local_norms_squared<kind>(matrix, uplo, norms);
        sqrt(norms);
        return norms;
      }
    else if constexpr(std::is_base_of_v<El::AbstractDistMatrix<El::BigFloat>,
                                        TMatrix>)
      {
        // TODO: in principle, each rank needs only the columns that it owns.
        // This could save some memory but introduce extra complexity.
        const auto comm = matrix.DistComm();
        std::vector<El::BigFloat> norms(norms_size<kind>(matrix), 0);
        add_local_norms_squared<kind>(matrix, uplo, norms);
        El::mpi::AllReduce(norms.data(), norms.size(), comm);
        sqrt(norms);
        return norms;
      }
    else
      {
        static_assert(always_false_v<TMatrix>,
                      "compute_norms(TMatrix) not implemented");
        return {};
      }
  }

  // Column norms for a Block_Matrix (vertically stacked blocks), where each columns
  // spans over all blocks stored on communicator comm.
  template <Matrix_Normalization_Kind kind>
  std::vector<El::BigFloat>
  compute_norms(const Block_Matrix &block_matrix, const El::mpi::Comm &comm)
  {
    static_assert(
      kind == Matrix_Normalization_Kind::COLUMNS,
      "compute_norms(Block_Matrix) is only valid for kind == COLUMNS");
    std::vector<El::BigFloat> norms(block_matrix.width, 0);
    for(const auto &matrix : block_matrix.blocks)
      {
        add_local_norms_squared<kind>(matrix, std::nullopt, norms);
      }
    El::mpi::AllReduce(norms.data(), norms.size(), comm);
    sqrt(norms);
    return norms;
  }
};

// Matrix_Normalizer implementations for different matrix classes

template <Matrix_Normalization_Kind kind, class TMatrix>
struct Abstract_Matrix_Normalizer
{
  using matrix_type = TMatrix;

protected:
  ~Abstract_Matrix_Normalizer() = default;

public:
  const int bits;
  const std::optional<El::UpperOrLower> uplo;
  const std::vector<El::BigFloat> norms;

  Abstract_Matrix_Normalizer(
    const matrix_type &matrix, const int bits,
    const std::optional<El::UpperOrLower> &uplo,
    const std::function<std::vector<El::BigFloat>(
      const matrix_type &, const std::optional<El::UpperOrLower> &)>
      compute_norms)
      : bits(bits), uplo(uplo), norms(compute_norms(matrix, uplo))
  {}

  // Normalize and multiply by 2^bits
  void normalize_and_shift(matrix_type &matrix,
                           const std::optional<El::UpperOrLower> &_uplo
                           = std::nullopt) const
  {
    // Parameter _uplo is added as a sanity check
    if(this->uplo != _uplo)
      {
        constexpr auto optional_to_string =
          [](
            const std::optional<El::UpperOrLower> &maybe_uplo) -> std::string {
          if(maybe_uplo.has_value())
            return std::to_string(El::UpperOrLowerToChar(maybe_uplo.value()));
          return "null";
        };
        LOGIC_ERROR("normalize_and_shift(): UpperOrLower value: ",
                    optional_to_string(_uplo),
                    " is different from stored in matrix normalizer: ",
                    optional_to_string(this->uplo));
      }
    Matrix_Normalizer_Util::update_local_elements(
      matrix,
      [&](const int iLoc, const int jLoc, El::BigFloat &element) {
        const auto &norm = get_local_element_norm(matrix, iLoc, jLoc);
        if(element == zero)
          return;
        // TODO is it actually possible due to rounding errors?
        ASSERT(norm != zero,
               "Vector with nonzero elements cannot have norm = 0",
               DEBUG_STRING(element));
        element /= norm;
        element <<= bits;
      },
      uplo);
  }
  // Reverse normalize_and_shift(), i.e divide by 2^bits and multiply by row/column norm.
  void restore(matrix_type &matrix) const
  {
    Matrix_Normalizer_Util::update_local_elements(
      matrix, [&](const int iLoc, const int jLoc, El::BigFloat &element) {
        const auto &norm = get_local_element_norm(matrix, iLoc, jLoc);
        element >>= bits;
        element *= norm;
      });
  }

private:
  const El::BigFloat &
  get_local_element_norm(matrix_type &matrix, int iLoc, int jLoc) const
  {
    const int row_or_col
      = Matrix_Normalizer_Util::global_row_or_col<kind>(matrix, iLoc, jLoc);
    return norms.at(row_or_col);
  }
  const El::BigFloat zero{0};
};

template <Matrix_Normalization_Kind kind>
struct Matrix_Normalizer
    : Abstract_Matrix_Normalizer<kind, El::Matrix<El::BigFloat>>
{
  Matrix_Normalizer(const El::Matrix<El::BigFloat> &matrix, int bits,
                    const std::optional<El::UpperOrLower> &uplo = std::nullopt)
      : Abstract_Matrix_Normalizer<kind, El::Matrix<El::BigFloat>>(
          matrix, bits, uplo,
          Matrix_Normalizer_Util::compute_norms<kind, El::Matrix<El::BigFloat>>)
  {}
};

template <Matrix_Normalization_Kind kind>
struct DistMatrix_Normalizer
    : Abstract_Matrix_Normalizer<kind, El::DistMatrix<El::BigFloat>>
{
  DistMatrix_Normalizer(const El::DistMatrix<El::BigFloat> &matrix,
                        const int bits,
                        const std::optional<El::UpperOrLower> &uplo
                        = std::nullopt)
      : Abstract_Matrix_Normalizer<kind, El::DistMatrix<El::BigFloat>>(
          matrix, bits, uplo,
          Matrix_Normalizer_Util::compute_norms<kind,
                                                El::DistMatrix<El::BigFloat>>)
  {}
};

template <Matrix_Normalization_Kind kind = Matrix_Normalization_Kind::COLUMNS>
struct Block_Matrix_Normalizer : Abstract_Matrix_Normalizer<kind, Block_Matrix>
{
  Block_Matrix_Normalizer(const Block_Matrix &matrix, const int bits,
                          const El::mpi::Comm &comm)
      : Abstract_Matrix_Normalizer<kind, Block_Matrix>(
          matrix, bits, std::nullopt,
          [&comm](const Block_Matrix &m, const auto &) {
            return Matrix_Normalizer_Util::compute_norms<kind>(m, comm);
          })
  {}
};

// Given a vector<TMatrix>, normalize each matrix independently.
template <class TMatrixNormalizer> struct Vector_TMatrix_Normalizer
{
  std::vector<TMatrixNormalizer> normalizers;

  template <class... TArgs>
  explicit Vector_TMatrix_Normalizer(
    const std::vector<typename TMatrixNormalizer::matrix_type> &matrices,
    const int bits, TArgs &&...args)
  {
    normalizers.reserve(matrices.size());
    for(const auto &matrix : matrices)
      normalizers.emplace_back(matrix, bits, std::forward<TArgs>(args)...);
  }

  // Normalize each block separately
  template <class... TArgs>
  void normalize_and_shift(
    std::vector<typename TMatrixNormalizer::matrix_type> &matrices,
    TArgs &&...args) const
  {
    ASSERT_EQUAL(matrices.size(), normalizers.size());
    for(size_t index = 0; index < matrices.size(); ++index)
      {
        normalizers.at(index).normalize_and_shift(
          matrices.at(index), std::forward<TArgs>(args)...);
      }
  }
};

// Normalize each block of vector<Block_Diagonal_Matrix>
template <Matrix_Normalization_Kind kind>
struct Block_Diagonal_Matrix_Normalizer
{
  using matrix_type = Block_Diagonal_Matrix;

  Vector_TMatrix_Normalizer<DistMatrix_Normalizer<kind>> normalizer;
  const std::vector<DistMatrix_Normalizer<kind>> &normalizers() const
  {
    return normalizer.normalizers;
  }

  Block_Diagonal_Matrix_Normalizer(
    const Block_Diagonal_Matrix &block_diagonal_matrix, int bits,
    const std::optional<El::UpperOrLower> &uplo = std::nullopt)
      : normalizer(block_diagonal_matrix.blocks, bits, uplo)
  {}

  // Normalize each block separately
  void normalize_and_shift(Block_Diagonal_Matrix &block_diagonal_matrix,
                           const std::optional<El::UpperOrLower> &uplo
                           = std::nullopt) const
  {
    normalizer.normalize_and_shift(block_diagonal_matrix.blocks, uplo);
  }
};

// Normalize each block of vector<Block_Diagonal_Matrix>
template <Matrix_Normalization_Kind kind>
using Vector_Block_Diagonal_Matrix_Normalizer
  = Vector_TMatrix_Normalizer<Block_Diagonal_Matrix_Normalizer<kind>>;

// Helper functions

// Helper to translate TMatrix -> TMatrix_Normalizer type
template <Matrix_Normalization_Kind kind, class TMatrix>
struct Matrix_Normalizer_Type
{
  static_assert(always_false_v<TMatrix>,
                "Matrix_Normalizer_Type for TMatrix is not specified.");
  using normalizer_type = void;
};
template <Matrix_Normalization_Kind kind>
struct Matrix_Normalizer_Type<kind, El::Matrix<El::BigFloat>>
{
  using normalizer_type = Matrix_Normalizer<kind>;
};
template <Matrix_Normalization_Kind kind>
struct Matrix_Normalizer_Type<kind, El::DistMatrix<El::BigFloat>>
{
  using normalizer_type = DistMatrix_Normalizer<kind>;
};
template <Matrix_Normalization_Kind kind>
struct Matrix_Normalizer_Type<kind, Block_Matrix>
{
  using normalizer_type = Block_Matrix_Normalizer<kind>;
};
template <Matrix_Normalization_Kind kind>
struct Matrix_Normalizer_Type<kind, Block_Diagonal_Matrix>
{
  using normalizer_type = Block_Diagonal_Matrix_Normalizer<kind>;
};
template <Matrix_Normalization_Kind kind>
struct Matrix_Normalizer_Type<kind, std::vector<Block_Diagonal_Matrix>>
{
  using normalizer_type = Vector_Block_Diagonal_Matrix_Normalizer<kind>;
};

template <Matrix_Normalization_Kind kind, class TMatrix, class... TArgs>
typename Matrix_Normalizer_Type<kind, TMatrix>::normalizer_type
normalize_and_shift(TMatrix &matrix, const int bits, const TArgs &...args)
{
  typename Matrix_Normalizer_Type<kind, TMatrix>::normalizer_type normalizer(
    matrix, bits, args...);
  normalizer.normalize_and_shift(matrix, args...);
  return normalizer;
}
// Specialize normalize_and_shift for Block_Matrix
// since Block_Matrix_Normalizer has an extra parameter in ctor.
// NB: ROWS are not allowed, so we specialize for COLUMNS only.
template <>
inline Block_Matrix_Normalizer<Matrix_Normalization_Kind::COLUMNS>
normalize_and_shift<Matrix_Normalization_Kind::COLUMNS>(
  Block_Matrix &matrix, const int bits, const El::mpi::Comm &comm)
{
  Block_Matrix_Normalizer<Matrix_Normalization_Kind::COLUMNS> normalizer(
    matrix, bits, comm);
  normalizer.normalize_and_shift(matrix);
  return normalizer;
}

// Restore the result of Q := P^T P
template <class TInputMatrix, class TOutputMatrix>
void restore_syrk_output(
  const El::UpperOrLower uplo,
  const Abstract_Matrix_Normalizer<Matrix_Normalization_Kind::COLUMNS,
                                   TInputMatrix> &syrk_input_normalizer,
  TOutputMatrix &syrk_output)
{
  namespace Util = Matrix_Normalizer_Util;
  const auto &column_norms = syrk_input_normalizer.norms;
  const auto bits = 2 * syrk_input_normalizer.bits;
  ASSERT_EQUAL(Util::global_height(syrk_output), column_norms.size());
  ASSERT_EQUAL(Util::global_width(syrk_output), column_norms.size());

  // Since each input vector (row or column) is normalized,
  // All normalized output elements should be in the range [-1-eps, 1+eps],
  // and all diagonal elements should be either 0 or (1 +- eps),
  // where the rounding error eps ~ O(2^{-input_bits}).
  // We use conservative upper bound eps = 2^{-input_bits/2}
  // for the assertion below.
  const El::BigFloat zero(0), one(1);
  const auto eps = one >> syrk_input_normalizer.bits / 2;
  const auto one_eps = one + eps;
  const auto minus_one_eps = -one_eps;

  Util::update_local_elements(
    syrk_output,
    [&](const int iLoc, const int jLoc, El::BigFloat &value) {
      const int i = Util::global_row(syrk_output, iLoc);
      const int j = Util::global_col(syrk_output, jLoc);
      value >>= bits;
      // Check normalization
      ASSERT(value >= minus_one_eps && value <= one_eps,
             "All normalized elements of Syrk output matrix "
             "should belong to the range [-1,1]. ",
             DEBUG_STRING(i), DEBUG_STRING(j), DEBUG_STRING(value),
             DEBUG_STRING(El::Abs(value) - one), DEBUG_STRING(eps));
      if(i == j)
        {
          ASSERT((column_norms.at(i) == zero && value == zero)
                   || El::Abs(value - one) <= eps,
                 "Normalized Syrk output matrix should have either zeros or "
                 "ones on diagonal. ",
                 DEBUG_STRING(i), DEBUG_STRING(j), DEBUG_STRING(value),
                 DEBUG_STRING(value - one), DEBUG_STRING(eps));
        }
      value *= column_norms.at(i) * column_norms.at(j);
    },
    uplo);
}

// Restore the result of Trmm call:
// X := op(L) * X | side = left
// X := X * op(L) | size = right
// L is (uplo) triangular
// op(L) == L or L^T (orientation)
template <Matrix_Normalization_Kind L_kind, Matrix_Normalization_Kind X_kind,
          class LMatrix, class XMatrix>
void restore_trmm_output(
  const El::LeftOrRight side, const El::UpperOrLower uplo,
  const El::Orientation orientation,
  const Abstract_Matrix_Normalizer<L_kind, LMatrix> &L_normalizer,
  const Abstract_Matrix_Normalizer<X_kind, XMatrix> &X_normalizer, XMatrix &X)
{
  namespace Util = Matrix_Normalizer_Util;

  ASSERT(L_normalizer.uplo.has_value(), "L matrix should be triangular");
  ASSERT_EQUAL(uplo, L_normalizer.uplo.value());
  ASSERT(!X_normalizer.uplo.has_value(),
         DEBUG_STRING(X_normalizer.uplo.value()),
         "Trmm requires one triangular and one rectangular matrix");

  const auto bits = L_normalizer.bits + X_normalizer.bits;

  auto left_kind = transpose(orientation, L_kind);
  auto right_kind = X_kind;

  const auto *left_row_norms = &L_normalizer.norms;
  const auto *right_column_norms = &X_normalizer.norms;
  if(side == El::RIGHT)
    {
      std::swap(left_row_norms, right_column_norms);
      std::swap(left_kind, right_kind);
    }
  ASSERT_EQUAL(left_kind, Matrix_Normalization_Kind::ROWS);
  ASSERT_EQUAL(right_kind, Matrix_Normalization_Kind::COLUMNS);

  ASSERT_EQUAL(X.Height(), left_row_norms->size());
  ASSERT_EQUAL(X.Width(), right_column_norms->size());

  // Since each input vector (row or column) is normalized,
  // All normalized output elements should be in the range [-1-eps, 1+eps],
  // where the rounding error eps ~ O(2^{-input_bits}).
  // We use conservative upper bound eps = 2^{-input_bits/2}
  // for the assertion below.
  const El::BigFloat one(1);
  const auto eps = one >> std::min(L_normalizer.bits, X_normalizer.bits) / 2;
  const auto one_eps = one + eps;
  const auto minus_one_eps = -one_eps;

  Util::update_local_elements(
    X, [&](const int iLoc, const int jLoc, El::BigFloat &value) {
      const int i = Util::global_row(X, iLoc);
      const int j = Util::global_col(X, jLoc);
      value >>= bits;
      // Check normalization
      ASSERT(value >= minus_one_eps && value <= one_eps,
             "All normalized elements of Trmm output matrix "
             "should belong to the range [-1,1]. ",
             DEBUG_STRING(i), DEBUG_STRING(j), DEBUG_STRING(value),
             DEBUG_STRING(El::Abs(value) - one), DEBUG_STRING(eps));
      value *= left_row_norms->at(i) * right_column_norms->at(j);
    });
}

// Trmm for (Block_Diagonal_Matrix L, std::vector<Block_Diagonal_Matrix> X)
template <Matrix_Normalization_Kind L_kind, Matrix_Normalization_Kind X_kind>
void restore_trmm_output(
  const El::LeftOrRight side, const El::UpperOrLower uplo,
  const El::Orientation orientation,
  const Block_Diagonal_Matrix_Normalizer<L_kind> &Ls_normalizer,
  const Vector_Block_Diagonal_Matrix_Normalizer<X_kind> &Xss_normalizer,
  std::vector<Block_Diagonal_Matrix> &Xss)
{
  const auto num_blocks = Ls_normalizer.normalizers().size();

  ASSERT_EQUAL(Xss_normalizer.normalizers.size(), Xss.size());
  for(size_t vec_index = 0; vec_index < Xss.size(); ++vec_index)
    {
      auto &Xs = Xss.at(vec_index);
      const auto &Xs_normalizer
        = Xss_normalizer.normalizers.at(vec_index).normalizer;
      ASSERT_EQUAL(Xs.blocks.size(), num_blocks, DEBUG_STRING(vec_index));
      ASSERT_EQUAL(Xs_normalizer.normalizers.size(), num_blocks,
                   DEBUG_STRING(vec_index));
      for(size_t block = 0; block < num_blocks; ++block)
        {
          const auto &L_normalizer = Ls_normalizer.normalizers().at(block);
          const auto &X_normalizer = Xs_normalizer.normalizers.at(block);
          auto &X = Xs.blocks.at(block);
          try
            {
              restore_trmm_output(side, uplo, orientation, L_normalizer,
                                  X_normalizer, X);
            }
          catch(std::exception &e)
            {
              RUNTIME_ERROR(DEBUG_STRING(vec_index), DEBUG_STRING(block),
                            e.what());
            }
        }
    }
}
