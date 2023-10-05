#pragma once

#include <flint/fmpz_mat.h>
#include <memory>
#include <El.hpp>

// RAII wrapper for FLINT's fmpz_mat_t (big integer matrix)
class Fmpz_Matrix
{
public:
  // NB: Do not call fmpz_mat_init or fmpz_mat_clear!
  // Memory is managed by automatically
  fmpz_mat_t fmpz_matrix{};

  Fmpz_Matrix();
  Fmpz_Matrix(int height, int width);

  Fmpz_Matrix(const Fmpz_Matrix &other);
  Fmpz_Matrix(Fmpz_Matrix &&other) noexcept;

  Fmpz_Matrix &operator=(const Fmpz_Matrix &other);
  Fmpz_Matrix &operator=(Fmpz_Matrix &&other) noexcept;

  virtual ~Fmpz_Matrix();

  // int - to match El::Matrix
  [[nodiscard]] int Height() const;
  [[nodiscard]] int Width() const;
  void Resize(int height, int width);

  [[nodiscard]] const fmpz *Get(int i, int j) const;
  void Set(int i, int j, const fmpz_t &value);
  fmpz *operator()(int i, int j);
  const fmpz *operator()(int i, int j) const;

  explicit Fmpz_Matrix(const El::Matrix<El::BigFloat> &input);
  explicit Fmpz_Matrix(const El::DistMatrix<El::BigFloat> &input);
  void ToBigFloatMatrix(El::Matrix<El::BigFloat> &output) const;
};
