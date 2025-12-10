#pragma once

#include "Abstract_Json_Element_Parser.hxx"

#include <filesystem>
#include <iostream>
#include <rapidjson/istreamwrapper.h>
#include <rapidjson/error/en.h>

template <class TParser>
void parse_json(std::istream &input_stream, TParser &parser,
                const std::filesystem::path &input_path = {})
{
  rapidjson::IStreamWrapper wrapper(input_stream);

  rapidjson::ParseResult res;
  try
    {
      rapidjson::Reader reader;
      res = reader.Parse(wrapper, parser);
    }
  catch(std::exception &e)
    {
      RUNTIME_ERROR("Failed to parse ", input_path,
                    ": offset=", wrapper.Tell(), ": ", e.what());
    }
  if(res.IsError())
    {
      RUNTIME_ERROR("Failed to parse ", input_path, ": offset=", res.Offset(),
                    ": error: ", rapidjson::GetParseError_En(res.Code()));
    }
}

template <class TParser>
void parse_json(const std::filesystem::path &input_path, TParser &parser)
{
  std::ifstream input_file(input_path);
  return parse_json(input_file, parser, input_path);
}
