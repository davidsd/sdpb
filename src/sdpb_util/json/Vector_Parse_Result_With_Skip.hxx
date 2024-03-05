#pragma once

#include <El.hpp>

// Result of parsing a vector,
// where some elements can be added to parsed_elements
// and others can be skipped
template <class TValue> struct Vector_Parse_Result_With_Skip
{
  std::vector<TValue> parsed_elements;
  // Indices of parsed elements
  std::vector<size_t> indices;
  // Total number of elements, including skipped
  size_t num_elements = 0;

  void reset()
  {
    parsed_elements.clear();
    indices.clear();
    num_elements = 0;
  }
  void add_element(TValue &&element, const size_t index)
  {
    ASSERT_EQUAL(index, num_elements);

    parsed_elements.emplace_back(std::forward<TValue>(element));
    indices.emplace_back(index);
    ++num_elements;
  }
  void skip_element(const size_t index)
  {
    ASSERT_EQUAL(index, num_elements);
    ++num_elements;
  }

  Vector_Parse_Result_With_Skip() = default;
  ~Vector_Parse_Result_With_Skip() = default;

  // Allow moving and prevent accidential copying

  Vector_Parse_Result_With_Skip(const Vector_Parse_Result_With_Skip &other)
    = delete;
  Vector_Parse_Result_With_Skip(Vector_Parse_Result_With_Skip &&other) noexcept
    = default;
  Vector_Parse_Result_With_Skip &
  operator=(const Vector_Parse_Result_With_Skip &other)
    = delete;
  Vector_Parse_Result_With_Skip &
  operator=(Vector_Parse_Result_With_Skip &&other) noexcept
    = default;
};