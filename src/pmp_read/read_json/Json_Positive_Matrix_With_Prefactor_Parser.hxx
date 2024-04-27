#pragma once

#include "Json_Damped_Rational_Parser.hxx"
#include "Json_Polynomial_Vector_Parser.hxx"
#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "sdpb_util/to_matrix.hxx"
#include "sdpb_util/json/Json_Vector_Parser_With_Skip.hxx"

using Vector_Of_Polynomial_Vectors = std::vector<Polynomial_Vector>;
using Matrix_Of_Polynomial_Vectors = std::vector<Vector_Of_Polynomial_Vectors>;

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
  std::optional<Damped_Rational> reduced_prefactor;
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
        polynomials_parser(skip,
                           [this](Matrix_Of_Polynomial_Vectors &&pvm) {
                             this->polynomials = std::move(pvm);
                           }),
        prefactor_parser(skip,
                         [this](Damped_Rational &&damped_rational) {
                           this->prefactor = std::move(damped_rational);
                         }),
        reduced_prefactor_parser(skip,
                                 [this](Damped_Rational &&damped_rational) {
                                   this->reduced_prefactor
                                     = std::move(damped_rational);
                                 }),
        sample_points_parser(
          skip,
          [this](std::vector<El::BigFloat> &&sample_points) {
            this->sample_points = std::move(sample_points);
          }),
        sample_scalings_parser(skip,
                               [this](std::vector<El::BigFloat> &&scalings) {
                                 this->sample_scalings = std::move(scalings);
                               }),
        reduced_sample_scalings_parser(
          skip,
          [this](std::vector<El::BigFloat> &&scalings) {
            this->reduced_sample_scalings = std::move(scalings);
          }),
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
  Json_Matrix_Of_Polynomial_Vectors_Parser polynomials_parser;
  Json_Damped_Rational_Parser prefactor_parser;
  Json_Damped_Rational_Parser reduced_prefactor_parser;
  Json_Float_Vector_Parser<El::BigFloat> sample_points_parser;
  Json_Float_Vector_Parser<El::BigFloat> sample_scalings_parser;
  Json_Float_Vector_Parser<El::BigFloat> reduced_sample_scalings_parser;
  Json_Polynomial_Vector_Parser bilinear_basis_parser;
  Json_Polynomial_Vector_Parser bilinear_basis_0_parser;
  Json_Polynomial_Vector_Parser bilinear_basis_1_parser;

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "polynomials")
      return polynomials_parser;
    if(key == "prefactor" || key == "DampedRational")
      return prefactor_parser;
    if(key == "reducedPrefactor")
      return reduced_prefactor_parser;
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
    const auto matrix = to_matrix(std::move(polynomials).value());
    // TODO add move ctor for Polynomial_Vector_Matrix?
    return Polynomial_Vector_Matrix(
      matrix, prefactor, reduced_prefactor, sample_points, sample_scalings,
      reduced_sample_scalings, bilinear_basis);
  }
  void clear_result() override
  {
    polynomials.reset();
    prefactor.reset();
    sample_points.reset();
    sample_scalings.reset();
    bilinear_basis.reset();
  }
  void reset_element_parsers(const bool skip) override
  {
    prefactor_parser.reset(skip);
    polynomials_parser.reset(skip);
    sample_points_parser.reset(skip);
    sample_scalings_parser.reset(skip);

    bilinear_basis_parser.reset(skip);
    bilinear_basis_0_parser.reset(skip);
    bilinear_basis_1_parser.reset(skip);
  }
};

using Json_Positive_Matrix_With_Prefactor_Array_Parser
  = Json_Vector_Parser_With_Skip<Json_Positive_Matrix_With_Prefactor_Parser>;