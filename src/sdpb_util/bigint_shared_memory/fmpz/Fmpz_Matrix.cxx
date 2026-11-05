#include "Fmpz_Matrix.hxx"
#include "fmpz_BigFloat_convert.hxx"

// Constructors

Fmpz_Matrix::Fmpz_Matrix() : Fmpz_Matrix(0, 0) {}

Fmpz_Matrix::Fmpz_Matrix(int height, int width)
{
  fmpz_mat_init(fmpz_matrix, height, width);
}
// Copy and move
Fmpz_Matrix::Fmpz_Matrix(const Fmpz_Matrix &other)
{
  if(Height() == other.Height() && Width() == other.Width())
    {
      fmpz_mat_set(fmpz_matrix, other.fmpz_matrix);
    }
  else
    {
      fmpz_mat_clear(fmpz_matrix);
      fmpz_mat_init_set(fmpz_matrix, other.fmpz_matrix);
    }
}
Fmpz_Matrix::Fmpz_Matrix(Fmpz_Matrix &&other) noexcept
{
  fmpz_mat_swap(fmpz_matrix, other.fmpz_matrix);
}
// Copy and move assignments
Fmpz_Matrix &Fmpz_Matrix::operator=(const Fmpz_Matrix &other)
{
  if(this == &other)
    return *this;
  fmpz_mat_clear(fmpz_matrix);
  fmpz_mat_init_set(fmpz_matrix, other.fmpz_matrix);
  return *this;
}

Fmpz_Matrix &Fmpz_Matrix::operator=(Fmpz_Matrix &&other) noexcept
{
  fmpz_mat_swap(fmpz_matrix, other.fmpz_matrix);
  return *this;
}

Fmpz_Matrix::~Fmpz_Matrix()
{
  fmpz_mat_clear(fmpz_matrix);
}

// Parameters

int Fmpz_Matrix::Height() const
{
  return fmpz_mat_nrows(fmpz_matrix);
}
int Fmpz_Matrix::Width() const
{
  return fmpz_mat_ncols(fmpz_matrix);
}
void Fmpz_Matrix::Resize(int height, int width)
{
  if(height == Height() && width == Width())
    return;
  fmpz_mat_clear(fmpz_matrix);
  fmpz_mat_init(fmpz_matrix, height, width);
}

// Member access

const fmpz *Fmpz_Matrix::Get(int i, int j) const
{
  return (*this)(i, j);
}
void Fmpz_Matrix::Set(int i, int j, const fmpz_t &value)
{
  fmpz_set((*this)(i, j), value);
}

fmpz *Fmpz_Matrix::operator()(int i, int j)
{
  return fmpz_mat_entry(fmpz_matrix, i, j);
}

const fmpz *Fmpz_Matrix::operator()(int i, int j) const
{
  return fmpz_mat_entry(fmpz_matrix, i, j);
}

// BigFloat convert

void Fmpz_Matrix::ToBigFloatMatrix(El::Matrix<El::BigFloat> &output) const
{
  int height = Height();
  int width = Width();
  output.Resize(height, width);

  for(int i = 0; i < height; ++i)
    for(int j = 0; j < width; ++j)
      {
        fmpz_t_to_BigFloat((*this)(i, j), output(i, j));
      }
}
Fmpz_Matrix::Fmpz_Matrix(const El::Matrix<El::BigFloat> &input)
    : Fmpz_Matrix(input.Height(), input.Width())
{
  int height = Height();
  int width = Width();
  for(int i = 0; i < height; ++i)
    for(int j = 0; j < width; ++j)
      {
        BigFloat_to_fmpz_t(input(i, j), (*this)(i, j));
      }
}
Fmpz_Matrix::Fmpz_Matrix(const El::DistMatrix<El::BigFloat> &input)
    : Fmpz_Matrix(input.Height(), input.Width())
{
  int height = Height();
  int width = Width();
  for(int i = 0; i < height; ++i)
    for(int j = 0; j < width; ++j)
      {
        BigFloat_to_fmpz_t(input.Get(i, j), (*this)(i, j));
      }
}
