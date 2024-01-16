#pragma once

#include "Json_Damped_Rational_Parser.hxx"
#include "Json_Polynomial_Parser.hxx"
#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "sdpb_util/to_matrix.hxx"
#include "sdpb_util/json/Json_Vector_Parser_With_Skip.hxx"

using Polynomial_Vector = std::vector<Polynomial>;
using Vector_Of_Polynomial_Vectors = std::vector<Polynomial_Vector>;
using Matrix_Of_Polynomial_Vectors = std::vector<Vector_Of_Polynomial_Vectors>;

using Json_Polynomial_Vector_Parser
  = Json_Vector_Parser<Json_Polynomial_Parser>;

using Json_Vector_Of_Polynomial_Vectors_Parser
  = Json_Vector_Parser<Json_Polynomial_Vector_Parser>;

// Matrix = Vector of Vectors <Polynomial_Vector>
using Json_Matrix_Of_Polynomial_Vectors_Parser
  = Json_Vector_Parser<Json_Vector_Of_Polynomial_Vectors_Parser>;

class Json_Positive_Matrix_With_Prefactor_Parser final
    : public Abstract_Json_Object_Parser<Polynomial_Vector_Matrix>
{
private:
  std::optional<Matrix_Of_Polynomial_Vectors> polynomials;
  std::optional<Damped_Rational> prefactor;

public:
  Json_Positive_Matrix_With_Prefactor_Parser(
    const bool skip,
    const std::function<void(Polynomial_Vector_Matrix &&)> &on_parsed,
    const std::function<void()> &on_skipped)
      : Abstract_Json_Object_Parser(skip, on_parsed, on_skipped),
        polynomials_parser(skip,
                           [this](Matrix_Of_Polynomial_Vectors &&pvm) {
                             polynomials = std::make_optional(std::move(pvm));
                           }),
        prefactor_parser(skip, [this](Damped_Rational &&prefactor) {
          this->prefactor = std::move(prefactor);
        })
  {}

private:
  Json_Matrix_Of_Polynomial_Vectors_Parser polynomials_parser;
  Json_Damped_Rational_Parser prefactor_parser;

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "polynomials")
      return polynomials_parser;
    if(key == "DampedRational")
      return prefactor_parser;
    El::RuntimeError(
      "Json_Positive_Matrix_With_Prefactor_Parser: unknown key=", key);
  }

public:
  value_type get_result() override
  {
    if(!polynomials.has_value())
      El::RuntimeError("rank=", El::mpi::Rank(), ": polynomials not found");
    const auto matrix = to_matrix(std::move(polynomials).value());
    // TODO add move ctor for Polynomial_Vector_Matrix?
    return Polynomial_Vector_Matrix(matrix, prefactor, std::nullopt,
                                    std::nullopt, std::nullopt);
  }
  void clear_result() override
  {
    polynomials.reset();
    prefactor.reset();
  }
  void reset_element_parsers(const bool skip) override
  {
    prefactor_parser.reset(skip);
    polynomials_parser.reset(skip);
  }
};

using Json_Positive_Matrix_With_Prefactor_Array_Parser
  = Json_Vector_Parser_With_Skip<Json_Positive_Matrix_With_Prefactor_Parser>;