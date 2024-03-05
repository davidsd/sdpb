#pragma once

#include <El.hpp>

#include <vector>

template <class T, class U>
El::Matrix<T> to_matrix(const std::vector<std::vector<U>> &elements,
                        std::function<T(const U &)> convert) noexcept(false)
{
  if(elements.empty())
    return {};

  El::Matrix<T> result(elements.size(), elements.at(0).size());
  for(int i = 0; i < result.Height(); ++i)
    {
      const auto &row = elements.at(i);
      ASSERT_EQUAL(row.size(), result.Width());
      for(int j = 0; j < result.Width(); ++j)
        {
          result.Set(i, j, convert(row.at(j)));
        }
    }

  return result;
}

template <class T>
El::Matrix<T>
to_matrix(const std::vector<std::vector<T>> &elements) noexcept(false)
{
  return to_matrix<T, T>(elements, [](const T &x) { return x; });
}