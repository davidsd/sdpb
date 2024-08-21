#pragma once

#include "Float.hxx"

namespace fs = std::filesystem;

// out.txt produced by SDPB
struct Parse_Sdpb_Out_Txt : boost::noncopyable
{
  std::string terminate_reason;
  std::map<std::string, Float> float_map;

  explicit Parse_Sdpb_Out_Txt(const fs::path &path)
  {
    CAPTURE(path);
    REQUIRE(is_regular_file(path));
    std::ifstream is(path);
    std::string line;
    while(std::getline(is, line))
      {
        // line format:
        // key = value;
        CAPTURE(line);
        std::vector<std::string> tokens;
        boost::split(tokens, line, boost::is_any_of("=;"));
        if(tokens.size() <= 1)
          break; // empty string (the last one)

        REQUIRE(tokens.size() == 3); // key, value, empty "" at the end
        auto key = tokens[0];
        auto value = tokens[1];
        boost::trim(key);
        boost::trim(value);
        if(key == "terminateReason")
          {
            // unquote string
            value = value.substr(1, value.size() - 2);
            terminate_reason = value;
          }
        else
          {
            float_map[key] = Float(value);
          }
      }
    REQUIRE(!terminate_reason.empty());
  }
};