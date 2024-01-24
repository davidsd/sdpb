#define RAPIDJSON_HAS_STDSTRING 1

#include <catch2/catch_amalgamated.hpp>

#include "pmp_read/read_json/Json_PMP_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/json/Json_Vector_Parser.hxx"
#include "test_util/diff.hxx"
#include "unit_tests/util/util.hxx"

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/error/en.h>

using Test_Util::REQUIRE_Equal::diff;

template <class TFloat = El::BigFloat>
using Json_Float_Vector_Parser_With_Skip
  = Json_Vector_Parser_With_Skip<Json_Float_Parser<TFloat>>;

// convert to json
namespace
{
  rapidjson::Value
  to_json(const El::BigFloat &value, rapidjson::Document &document)
  {
    std::stringstream ss;
    set_stream_precision(ss);
    ss.str("");
    ss << value;
    rapidjson::Value result(rapidjson::kStringType);
    result.SetString(ss.str(), document.GetAllocator());
    return result;
  }

  template <class T>
  rapidjson::Value
  to_json(const std::vector<T> &values, rapidjson::Document &document)
  {
    rapidjson::Value result(rapidjson::kArrayType);
    auto &arr = result.SetArray();
    for(const auto &val : values)
      {
        arr.PushBack(to_json(val, document), document.GetAllocator());
      }
    return result;
  }

  template <class T> std::string to_json_string(const std::vector<T> &values)
  {
    rapidjson::Document document;
    document.SetArray() = to_json(values, document).GetArray();
    rapidjson::StringBuffer buffer;
    rapidjson::Writer writer(buffer);
    document.Accept(writer);
    return buffer.GetString();
  }
  std::string to_json_string(const Polynomial &poly)
  {
    return to_json_string(poly.coefficients);
  }
}

// parse json
namespace
{
  template <class TParser>
  void json_parse(const std::string &json_string, TParser &parser)
  {
    std::stringstream ss(json_string);
    rapidjson::IStreamWrapper wrapper(ss);
    rapidjson::ParseResult res;
    try
      {
        rapidjson::Reader reader;
        res = reader.Parse(wrapper, parser);
      }
    catch(std::exception &e)
      {
        CAPTURE(wrapper.Tell());
        FAIL(e.what());
        throw;
      }
    if(res.IsError())
      {
        CAPTURE(res.Offset());
        CAPTURE(rapidjson::GetParseError_En(res.Code()));
        CAPTURE(std::string(json_string).substr(res.Offset(), 10));
        FAIL("Parser failed");
      }
  }

  template <class TFloat = El::BigFloat>
  std::vector<TFloat> json_to_vector(const std::string &json_string)
  {
    std::vector<TFloat> output_values;

    bool skip = false;
    auto on_parsed = [&output_values](auto &&parse_result) {
      output_values = std::forward<std::vector<El::BigFloat>>(parse_result);
    };
    Json_Float_Vector_Parser<El::BigFloat> parser(skip, on_parsed);

    json_parse(json_string, parser);
    return output_values;
  }

  template <class TFloat = El::BigFloat>
  std::vector<TFloat> json_to_vector_with_skip(
    const std::string &json_string,
    const std::function<bool(size_t index)> &skip_element)
  {
    std::vector<TFloat> output_values;

    bool skip = false;
    auto on_parsed = [&output_values](auto &&parse_result) {
      output_values = std::move(parse_result.parsed_elements);
    };
    auto on_skipped = [] {};
    Json_Float_Vector_Parser_With_Skip<TFloat> parser(
      skip, on_parsed, on_skipped, [&skip_element](size_t index) {
        const bool skip_index = skip_element(index);
        return skip_index;
      });

    json_parse(json_string, parser);
    return output_values;
  }

  Polynomial json_to_Polynomial(const std::string &json_string)
  {
    Polynomial output;
    output.coefficients.clear();

    auto on_parsed = [&output](auto &&parse_result) {
      output = std::forward<Polynomial>(parse_result);
    };
    const bool skip = false;

    Json_Polynomial_Parser parser(skip, on_parsed);

    json_parse(json_string, parser);
    return output;
  }

  template <class TFloat = El::BigFloat>
  std::vector<std::vector<TFloat>> json_to_vector_vector_with_skip(
    const std::string &json_string,
    const std::function<bool(size_t index)> &outer_skip_func,
    const std::function<bool(size_t index)> &inner_skip_func)
  {
    // Parse JSON
    std::vector<std::vector<El::BigFloat>> output_values_with_skip;

    auto parser = Json_Vector_Parser_With_Skip<
      Json_Float_Vector_Parser_With_Skip<El::BigFloat>>(
      false,
      [&output_values_with_skip](auto &&result) {
        for(auto &inner_result : result.parsed_elements)
          {
            output_values_with_skip.emplace_back(inner_result.parsed_elements);
          }
      },
      [] {}, outer_skip_func, inner_skip_func);

    json_parse(json_string, parser);
    return output_values_with_skip;
  }
}

TEST_CASE("json")
{
  if(El::mpi::Rank() != 0)
    return;

  SECTION("BigFloat_Vector_Parser")
  {
    const std::vector<El::BigFloat> input_values{0.1, -2, 34.5};
    CAPTURE(input_values);

    const auto json_string = to_json_string(input_values);
    CAPTURE(json_string);
    std::vector<El::BigFloat> output_values
      = json_to_vector<El::BigFloat>(json_string);

    DIFF(input_values, output_values);
  }

  SECTION("Polynomial")
  {
    Polynomial input_poly;
    input_poly.coefficients = {0.1, -2, 34.5};
    const auto json_string = to_json_string(input_poly);
    CAPTURE(json_string);
    Polynomial output_poly = json_to_Polynomial(json_string);

    DIFF(input_poly, output_poly);
  }

  SECTION("BigFloat_Vector_Parser_With_Skip")
  {
    const std::vector<El::BigFloat> input_values{0.1, -2, 34.5};
    std::vector<El::BigFloat> input_even;
    std::vector<El::BigFloat> input_odd;
    for(size_t i = 0; i < input_values.size(); ++i)
      {
        if(i % 2 == 0)
          input_even.push_back(input_values.at(i));
        else
          input_odd.push_back(input_values.at(i));
      }

    const auto json_string = to_json_string(input_values);
    CAPTURE(json_string);

    {
      std::vector<El::BigFloat> output_even = json_to_vector_with_skip(
        json_string, [](size_t index) { return index % 2 != 0; });
      CAPTURE(input_even);
      CAPTURE(output_even);
      DIFF(input_even, output_even);
    }
    {
      std::vector<El::BigFloat> output_odd = json_to_vector_with_skip(
        json_string, [](size_t index) { return index % 2 == 0; });
      CAPTURE(input_odd);
      CAPTURE(output_odd);
      DIFF(input_odd, output_odd);
    }
  }

  SECTION("Json_Vector_Vector_With_Skip")
  {
    INFO("Check parsing vector of vectors with element skipping "
         "in both outer and inner vectors.");
    // Random vector of vectors
    std::vector<std::vector<El::BigFloat>> input_values;
    for(size_t i = 0; i < 5; ++i)
      {
        input_values.push_back(Test_Util::random_vector(5 + i));
      }

    CAPTURE(input_values);

    auto json_string = to_json_string(input_values);

    // Skip odd or even elements ofouter and inner vector

    for(const bool outer_skip_even : {false, true})
      {
        auto outer_skip_func = [outer_skip_even](size_t index) {
          return (index % 2 == 0) == outer_skip_even;
        };
        CAPTURE(outer_skip_even);
        for(const bool inner_skip_even : {false, true})
          {
            CAPTURE(inner_skip_even);
            auto inner_skip_func = [inner_skip_even](size_t index) {
              return (index % 2 == 0) == inner_skip_even;
            };

            // Select values from input
            std::vector<std::vector<El::BigFloat>> input_values_with_skip;
            for(size_t i = 0; i < input_values.size(); ++i)
              {
                if(outer_skip_func(i))
                  continue;
                auto &inner_vec = input_values_with_skip.emplace_back();
                for(size_t j = 0; j < input_values.at(i).size(); ++j)
                  {
                    if(inner_skip_func(j))
                      continue;
                    inner_vec.push_back(input_values.at(i).at(j));
                  }
              }
            CAPTURE(input_values_with_skip);
            REQUIRE(!input_values_with_skip.empty());

            // Parse JSON
            std::vector<std::vector<El::BigFloat>> output_values_with_skip
              = json_to_vector_vector_with_skip(json_string, outer_skip_func,
                                                inner_skip_func);

            // Check result
            CAPTURE(output_values_with_skip);
            DIFF(input_values_with_skip, output_values_with_skip);
          }
      }
  }
}
