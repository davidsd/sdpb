#include "Json_PMP_Parser.hxx"
#include "sdpb_util/json/Vector_Parse_Result_With_Skip.hxx"

Json_PMP_Parser::Json_PMP_Parser(
  const std::function<bool(size_t matrix_index)> &should_parse_matrix,
  const std::function<void(value_type &&result)> &on_parsed)
    : Abstract_Json_Object_Parser(false, on_parsed, [] {}),
      objective_parser(
        // don't skip objective:
        false,
        [this](std::vector<El::BigFloat> &&result) {
          this->result.objective = std::move(result);
        },
        [] { LOGIC_ERROR(R"(Skipping "objective" not allowed)"); }),
      normalization_parser(
        // don't skip normalization:
        false,
        // accept normalization vector:
        [this](std::vector<El::BigFloat> &&result) {
          this->result.normalization = std::move(result);
        },
        [] { LOGIC_ERROR(R"(Skipping "normalization" not allowed)"); }),
      matrices_parser(
        // don't skip matrices array:
        false,
        // accept matrix
        [this](
          Vector_Parse_Result_With_Skip<Polynomial_Vector_Matrix> &&result) {
          this->result.num_matrices = result.num_elements;
          assert(result.indices.size() == result.parsed_elements.size());
          for(size_t i = 0; i < result.indices.size(); ++i)
            {
              this->result.parsed_matrices.emplace(
                result.indices.at(i), std::move(result.parsed_elements.at(i)));
            }
        },
        [] {
          LOGIC_ERROR(
            R"(Skipping "PositiveMatrixWithPrefactorArray" not allowed)");
        },
        // Skip some matrices according to their indices:
        [&should_parse_matrix](size_t index) {
          return !should_parse_matrix(index);
        })
{}
Abstract_Json_Reader_Handler &
Json_PMP_Parser::element_parser(const std::string &key)
{
  if(key == "objective")
    return objective_parser;
  if(key == "normalization")
    return normalization_parser;
  if(key == "PositiveMatrixWithPrefactorArray")
    return matrices_parser;
  RUNTIME_ERROR("Json_PMP_Parser: Unexpected key=", key);
}
void Json_PMP_Parser::clear_result()
{
  result.objective.clear();
  result.normalization.clear();
  result.num_matrices = 0;
  result.parsed_matrices.clear();
  result.path.clear();
}
Json_PMP_Parser::value_type Json_PMP_Parser::get_result()
{
  return std::move(result);
}
void Json_PMP_Parser::reset_element_parsers(const bool skip)
{
  objective_parser.reset(skip);
  normalization_parser.reset(skip);
  matrices_parser.reset(skip);
}
