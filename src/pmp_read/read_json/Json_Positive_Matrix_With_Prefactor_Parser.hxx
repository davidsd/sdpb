#pragma once

#include "Json_Damped_Rational_Parser.hxx"
#include "Json_Polynomial_Power_Product_Parser.hxx"
#include "Json_Polynomial_Parser.hxx"
#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "sdpb_util/to_matrix.hxx"

using Vector_Of_Polynomial_Vectors = std::vector<Polynomial_Vector>;
using Matrix_Of_Polynomial_Vectors = std::vector<Vector_Of_Polynomial_Vectors>;

class Json_Positive_Matrix_With_Prefactor_Parser final
    : public Abstract_Json_Object_Parser<Polynomial_Vector_Matrix>
{
private:
  std::optional<Matrix_Of_Polynomial_Vectors> polynomials;
  std::optional<Damped_Rational> prefactor;
  std::optional<Damped_Rational> reduced_prefactor;
  std::optional<std::vector<Polynomial_Power_Product>> preconditioning_vector;
  std::optional<std::vector<El::BigFloat>> sample_points;
  std::optional<std::vector<El::BigFloat>> sample_scalings;
  std::optional<std::vector<El::BigFloat>> reduced_sample_scalings;
  std::optional<std::array<Polynomial_Vector, 2>> bilinear_basis;

public:
  Json_Positive_Matrix_With_Prefactor_Parser(
    const bool skip,
    const std::function<void(Polynomial_Vector_Matrix &&)> &on_parsed,
    const std::function<void()> &on_skipped)
      : Abstract_Json_Object_Parser(skip, on_parsed, on_skipped),

#define ELEMENT_PARSER_CTOR(element_name)                                     \
  element_name##_parser(                                                      \
    skip, [this](auto &&value) { this->element_name = std::move(value); })

        ELEMENT_PARSER_CTOR(polynomials),
        ELEMENT_PARSER_CTOR(prefactor),
        ELEMENT_PARSER_CTOR(reduced_prefactor),
        ELEMENT_PARSER_CTOR(preconditioning_vector),
        ELEMENT_PARSER_CTOR(sample_points),
        ELEMENT_PARSER_CTOR(sample_scalings),
        ELEMENT_PARSER_CTOR(reduced_sample_scalings),

#undef ELEMENT_PARSER_CTOR

        bilinear_basis_parser(skip,
                              [this](Polynomial_Vector &&value) {
                                if(!this->bilinear_basis.has_value())
                                  this->bilinear_basis = {value, value};
                              }),
        bilinear_basis_0_parser(skip,
                                [this](Polynomial_Vector &&value) {
                                  if(!this->bilinear_basis.has_value())
                                    this->bilinear_basis.emplace();
                                  this->bilinear_basis.value()[0]
                                    = std::move(value);
                                }),
        bilinear_basis_1_parser(skip, [this](Polynomial_Vector &&value) {
          if(!this->bilinear_basis.has_value())
            this->bilinear_basis.emplace();
          this->bilinear_basis.value()[1] = std::move(value);
        })
  {}

private:
  using BigFloat_Vector_Parser = Json_Vector_Parser<Json_BigFloat_Parser>;
  using Damped_Rational_Parser
    = Json_Damped_Rational_Parser;
  using Polynomial_Power_Product_Vector_Parser=  Json_Vector_Parser<Json_Polynomial_Power_Product_Parser>;
  using Polynomial_Vector_Parser
    = Json_Vector_Parser<Json_Polynomial_Parser>;
  using Vector_Of_Polynomial_Vectors_Parser
    = Json_Vector_Parser<Polynomial_Vector_Parser>;
  // Matrix = Vector of Vectors <Polynomial_Vector>
  using Matrix_Of_Polynomial_Vectors_Parser
    = Json_Vector_Parser<Vector_Of_Polynomial_Vectors_Parser>;

  Matrix_Of_Polynomial_Vectors_Parser polynomials_parser;
  Damped_Rational_Parser prefactor_parser;
  Damped_Rational_Parser reduced_prefactor_parser;
  Polynomial_Power_Product_Vector_Parser preconditioning_vector_parser;
  BigFloat_Vector_Parser sample_points_parser;
  BigFloat_Vector_Parser sample_scalings_parser;
  BigFloat_Vector_Parser reduced_sample_scalings_parser;
  Polynomial_Vector_Parser bilinear_basis_parser;
  Polynomial_Vector_Parser bilinear_basis_0_parser;
  Polynomial_Vector_Parser bilinear_basis_1_parser;

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "polynomials")
      return polynomials_parser;
    if(key == "prefactor" || key == "DampedRational")
      return prefactor_parser;
    if(key == "reducedPrefactor")
      return reduced_prefactor_parser;
    if(key == "preconditioningVector")
      return preconditioning_vector_parser;

    if(key == "samplePoints")
      return sample_points_parser;
    if(key == "sampleScalings")
      return sample_scalings_parser;
    if(key == "reducedSampleScalings")
      return reduced_sample_scalings_parser;

    if(key == "bilinearBasis")
      return bilinear_basis_parser;
    if(key == "bilinearBasis_0")
      return bilinear_basis_0_parser;
    if(key == "bilinearBasis_1")
      return bilinear_basis_1_parser;
    RUNTIME_ERROR("unknown key=", key);
  }

public:
  value_type get_result() override
  {
    ASSERT(polynomials.has_value(), "polynomials not found");
    const Simple_Matrix matrix(std::move(polynomials).value());
    // TODO add move ctor for Polynomial_Vector_Matrix?
    return Polynomial_Vector_Matrix(
      matrix, prefactor, reduced_prefactor, preconditioning_vector,
      sample_points, sample_scalings, reduced_sample_scalings, bilinear_basis);
  }
  void clear_result() override
  {
    polynomials.reset();
    prefactor.reset();
    reduced_prefactor.reset();
    preconditioning_vector.reset();

    sample_points.reset();
    sample_scalings.reset();
    reduced_sample_scalings.reset();
    bilinear_basis.reset();
  }
  void reset_element_parsers(const bool skip) override
  {
    polynomials_parser.reset(skip);
    prefactor_parser.reset(skip);
    reduced_prefactor_parser.reset(skip);
    preconditioning_vector_parser.reset(skip);

    sample_points_parser.reset(skip);
    sample_scalings_parser.reset(skip);
    reduced_sample_scalings_parser.reset(skip);

    bilinear_basis_parser.reset(skip);
    bilinear_basis_0_parser.reset(skip);
    bilinear_basis_1_parser.reset(skip);
  }
};
