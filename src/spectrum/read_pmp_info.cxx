#include "pmp/PMP_Info.hxx"
#include "pmp_read/read_json/Json_Polynomial_Power_Product_Parser.hxx"
#include "sdpb_util/Archive_Reader.hxx"
#include "sdpb_util/Timers/Timers.hxx"
#include "sdpb_util/json/Abstract_Json_Object_Parser.hxx"
#include "sdpb_util/json/Json_Damped_Rational_Parser.hxx"
#include "sdpb_util/json/Json_Float_Parser.hxx"
#include "sdpb_util/json/Json_Skip_Element_Parser.hxx"
#include "sdpb_util/json/Json_String_Parser.hxx"
#include "sdpb_util/json/Json_UInt64_Parser.hxx"
#include "sdpb_util/json/Json_Vector_Parser_With_Skip.hxx"
#include "sdpb_util/json/parse_json.hxx"

namespace fs = std::filesystem;

class Json_PVM_Info_Parser final : public Abstract_Json_Object_Parser<PVM_Info>
{
  // TODO store separate fields an assert that they are all initialized?
  PVM_Info result;

  Json_UInt64_Parser block_index_parser;
  Json_String_Parser block_path_parser;
  Json_UInt64_Parser dim_parser;
  Json_Damped_Rational_Parser prefactor_parser;
  Json_Damped_Rational_Parser reduced_prefactor_parser;
  Json_Vector_Parser<Json_Polynomial_Power_Product_Parser>
    preconditioning_vector_parser;
  Json_Vector_Parser<Json_BigFloat_Parser> sample_points_parser;
  Json_Vector_Parser<Json_BigFloat_Parser> sample_scalings_parser;
  Json_Vector_Parser<Json_BigFloat_Parser> reduced_sample_scalings_parser;
  Json_Skip_Element_Parser skip_element_parser{};

public:
  Json_PVM_Info_Parser(const bool skip,
                       const std::function<void(PVM_Info &&)> &on_parsed,
                       const std::function<void()> &on_skipped)
      : Abstract_Json_Object_Parser(skip, on_parsed, on_skipped),

#define ELEMENT_PARSER_CTOR(element_name)                                     \
  element_name##_parser(                                                      \
    skip, [this](auto &&value) { this->result.element_name = value; })
        ELEMENT_PARSER_CTOR(block_index),
        ELEMENT_PARSER_CTOR(block_path),
        ELEMENT_PARSER_CTOR(dim),
        ELEMENT_PARSER_CTOR(prefactor),
        ELEMENT_PARSER_CTOR(reduced_prefactor),
        ELEMENT_PARSER_CTOR(preconditioning_vector),
        ELEMENT_PARSER_CTOR(sample_points),
        ELEMENT_PARSER_CTOR(sample_scalings),
        ELEMENT_PARSER_CTOR(reduced_sample_scalings)
#undef ELEMENT_PARSER_CTOR
  {}

public:
  value_type get_result() override { return result; }

protected:
  Abstract_Json_Reader_Handler &element_parser(const std::string &key) override
  {
    if(key == "index")
      return block_index_parser;
    if(key == "path")
      return block_path_parser;

    if(key == "dim")
      return dim_parser;

    if(key == "prefactor")
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

    PRINT_WARNING("unknown pmp_info JSON key=\"", key, "\"");
    return skip_element_parser;
  }

public:
  void reset_element_parsers(const bool skip) override
  {
    block_index_parser.reset(skip);
    block_path_parser.reset(skip);
    dim_parser.reset(skip);
    prefactor_parser.reset(skip);
    reduced_prefactor_parser.reset(skip);
    sample_points_parser.reset(skip);
    sample_scalings_parser.reset(skip);
    reduced_sample_scalings_parser.reset(skip);
    preconditioning_vector_parser.reset(skip);
  }
  void clear_result() override { result = PVM_Info(); }
};

PMP_Info read_pmp_info(const std::filesystem::path &input_path, Timers &timers)
{
  Scoped_Timer timer(timers, "read_pmp_info");

  Vector_Parse_Result_With_Skip<PVM_Info> parse_result;

  const bool skip = false;
  const auto on_parsed = [&parse_result](auto &&value) {
    parse_result = std::forward<decltype(value)>(value);
  };
  const auto on_skipped = [] {};

  // Round-robin block distribution among ranks.
  // TODO: use something smarter, based on block sizes?
  // e.g. rank=0 reads block sizes from pmp_info.json,
  // decides on block distribution,
  // and then each rank reads its blocks pmp_info.json again.
  // Note also that we are reading a single file with multiple processes.
  // Maybe it's better to read it on rank=0 and send data to other ranks?
  const auto should_skip_block = [](const size_t index) {
    return index % El::mpi::Size() != El::mpi::Rank();
  };

  // pmp_info.json contains an array of PVM_Info objects
  Json_Vector_Parser_With_Skip<Json_PVM_Info_Parser> parser(
    skip, on_parsed, on_skipped, should_skip_block);

  // read pmp_info.json from regular file
  if(fs::exists(input_path))
    {
      parse_json(input_path, parser);
    }
  // read pmp_info.json from sdp.zip archive
  else
    {
      const auto parent_path = input_path.parent_path();
      if(fs::exists(parent_path) && !fs::is_directory(parent_path))
        {
          const auto filename = input_path.filename();
          try
            {
              Archive_Reader reader(parent_path);
              bool found_pmp_info = false;
              while(reader.next_entry())
                {
                  if(filename == archive_entry_pathname(reader.entry_ptr))
                    {
                      std::istream stream(&reader);
                      parse_json(stream, parser, input_path);
                      found_pmp_info = true;
                      break;
                    }
                }
              ASSERT(found_pmp_info, "Unable to find ", filename, " in ",
                     parent_path);
            }
          catch(const std::exception &e)
            {
              RUNTIME_ERROR("Failed to read ", filename, " from archive ",
                            parent_path, ": ", e.what());
            }
        }
      else
        {
          RUNTIME_ERROR("Cannot find file: ", input_path);
        }
    }

  // Validate data
  for(size_t i = 0; i < parse_result.parsed_elements.size(); ++i)
    {
      const auto &index = parse_result.indices.at(i);
      auto &block = parse_result.parsed_elements.at(i);
      // if "index" was missing from JSON
      if(block.block_index < 0)
        block.block_index = index;
      ASSERT(block.block_index == index, "For block_", index,
             ", read \"index\"=", block.block_index, " expected ", index);
      block.validate(parse_result.num_elements);
    }

  return PMP_Info(parse_result.num_elements, parse_result.parsed_elements);
}
