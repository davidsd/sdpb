#pragma once

#include "Json_Vector_Parser.hxx"
#include "sdpb_util/to_matrix.hxx"

#include <El.hpp>

#include <vector>

// Parses JSON array of arrays into El::Matrix<T> (row by row),
// each element T is parsed with TElementParser
template <class TMatrixElementParser>
class Json_Matrix_Parser final
    : public Abstract_Json_Vector_Parser<
        El::Matrix<typename TMatrixElementParser::value_type>,
        Json_Vector_Parser<TMatrixElementParser>>
{
public:
  virtual ~Json_Matrix_Parser() = default;
  using base_type = Abstract_Json_Vector_Parser<
    El::Matrix<typename TMatrixElementParser::value_type>,
    Json_Vector_Parser<TMatrixElementParser>>;
  // std::vector<T>
  using element_type = typename base_type::element_type;
  // El::Matrix<T>
  using value_type = El::Matrix<typename TMatrixElementParser::value_type>;

private:
  // std::vector<std::vector<T>>, to be converted to El::Matrix<T>
  std::vector<typename base_type::element_type> elements;

public:
  Json_Matrix_Parser(
    const bool skip, const std::function<void(value_type &&)> &on_parsed,
    const std::function<void()> &on_skipped = [] {})
      : base_type(skip, on_parsed, on_skipped)
  {}

  // Parsed single row
  void on_element_parsed(element_type &&value, size_t index) override
  {
    ASSERT(index == elements.size(), "index=", index, ", expected ",
           elements.size());
    elements.push_back(std::move(value));
  }

private:
  void clear_result() override { elements.clear(); }
  value_type get_result() override { return to_matrix(elements); }
};
