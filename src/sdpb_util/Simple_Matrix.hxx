#pragma once

#include "assert.hxx"

#include <vector>

// Simple interface for a fixed-size matrix
template <class T> class Simple_Matrix
{
private:
  const size_t height;
  const size_t width;
  std::vector<T> data;

public:
  // Constructors

  Simple_Matrix(const size_t height, const size_t width)
      : height(height), width(width), data(height * width)
  {}

  Simple_Matrix() : Simple_Matrix(0, 0) {}

  explicit Simple_Matrix(const std::vector<std::vector<T>> &rows)
      : Simple_Matrix(rows.size(), rows.empty() ? 0 : rows[0].size())
  {
    for(size_t i = 0; i < height; i++)
      {
        const auto &row = rows[i];
        ASSERT_EQUAL(row.size(), width, DEBUG_STRING(i),
                     "Input has different row lengths!");
        for(size_t j = 0; j < width; j++)
          {
            (*this)(i, j) = row[j];
          }
      }
  }

  // El::Matrix-like interface

  [[nodiscard]] size_t Height() const { return height; }

  [[nodiscard]] size_t Width() const { return width; }

  T &operator()(const size_t i, const size_t j)
  {
    ASSERT(i < height, DEBUG_STRING(i), DEBUG_STRING(height));
    ASSERT(j < width, DEBUG_STRING(j), DEBUG_STRING(width));
    return data.at(i + j * height);
  }

  const T &operator()(const size_t i, const size_t j) const
  {
    ASSERT(i < height, DEBUG_STRING(i), DEBUG_STRING(height));
    ASSERT(j < width, DEBUG_STRING(j), DEBUG_STRING(width));
    return data.at(i + j * height);
  }
};