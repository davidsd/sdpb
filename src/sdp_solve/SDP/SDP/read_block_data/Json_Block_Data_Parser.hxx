#pragma once

#include "Block_Data_Parse_Result.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/json/Json_Matrix_Parser.hxx"

using namespace std::string_literals;
struct Json_Block_Data_Parser final
    : Abstract_Json_Object_Parser<Block_Data_Parse_Result>
{
private:
  Block_Data_Parse_Result result;

  Json_Float_Vector_Parser<El::BigFloat> c_parser;
  Json_Matrix_Parser<Json_Float_Parser<El::BigFloat>> B_parser;
  Json_Matrix_Parser<Json_Float_Parser<El::BigFloat>>
    bilinear_bases_even_parser;
  Json_Matrix_Parser<Json_Float_Parser<El::BigFloat>> bilinear_bases_odd_parser;

  // Abstract_Json_Object_Parser implementation
public:
  Block_Data_Parse_Result get_result() override { return std::move(result); }

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "c")
      return c_parser;
    if(key == "B")
      return B_parser;
    if(key == "bilinear_bases_even")
      return bilinear_bases_even_parser;
    if(key == "bilinear_bases_odd")
      return bilinear_bases_odd_parser;
    RUNTIME_ERROR("Unknown key: '", key, "'");
  }

public:
  void reset_element_parsers(bool skip) override
  {
    c_parser.reset(skip);
    B_parser.reset(skip);
    bilinear_bases_even_parser.reset(skip);
    bilinear_bases_odd_parser.reset(skip);
  }
  void clear_result() override { result.clear(); }

public:
  explicit Json_Block_Data_Parser(
    const std::function<void(Block_Data_Parse_Result &&)> &on_parsed)
      : Abstract_Json_Object_Parser(false, on_parsed, [] {}),
        c_parser(false,
                 [this](std::vector<El::BigFloat> &&c) {
                   this->result.c = std::move(c);
                 }),
        B_parser(false,
                 [this](El::Matrix<El::BigFloat> &&B) {
                   this->result.B = std::move(B);
                 }),
        bilinear_bases_even_parser(
          false,
          [this](El::Matrix<El::BigFloat> &&bilinear_bases_even) {
            this->result.bilinear_bases_even = std::move(bilinear_bases_even);
          }),
        bilinear_bases_odd_parser(
          false, [this](El::Matrix<El::BigFloat> &&bilinear_bases_odd) {
            this->result.bilinear_bases_odd = std::move(bilinear_bases_odd);

            // See #124 sdpb fails to parse block_data with empty bilinear_bases_odd
            if(this->result.bilinear_bases_odd.Width() == 0)
              {
                this->result.bilinear_bases_odd.Resize(
                  result.bilinear_bases_odd.Height(), 1);
              }
          })
  {}
};
