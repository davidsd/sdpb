#pragma once

#include "sdpb_util/assert.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/json/Json_Vector_Parser.hxx"

using namespace std::string_literals;
struct Json_Points_Parser final
    : Abstract_Json_Object_Parser<std::vector<std::vector<El::BigFloat>>>
{
private:
  std::vector<std::vector<El::BigFloat>> result;
  Json_Vector_Parser<Json_Float_Vector_Parser<El::BigFloat>> points_parser;

  // Abstract_Json_Object_Parser implementation
public:
  std::vector<std::vector<El::BigFloat>> get_result() override
  {
    return std::move(result);
  }

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "points")
      return points_parser;
    RUNTIME_ERROR("Unknown key: '", key, "'");
  }

public:
  void reset_element_parsers(bool skip) override { points_parser.reset(skip); }
  void clear_result() override { result.clear(); }

  explicit
  Json_Points_Parser(const std::function<void(value_type &&result)> &on_parsed)
      : Abstract_Json_Object_Parser(false, on_parsed, [] {}),
        points_parser(false, [this](value_type &&points) {
          this->result = std::move(points);
        })
  {}
};
