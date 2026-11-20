#include "read_normalization.hxx"

#include "sdpb_util/Archive_Reader.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <rapidjson/document.h>
#include <rapidjson/istreamwrapper.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <optional>

namespace fs = std::filesystem;

namespace
{
  // Each value keeps all digits of normalization.json, which can be more
  // than SDPB's precision: z.txt needs them (sdpb/z_normalization.hxx).
  std::vector<Boost_Float>
  read_normalization_stream(std::istream &normalization_stream)
  {
    rapidjson::IStreamWrapper wrapper(normalization_stream);
    rapidjson::Document d;
    d.ParseStream(wrapper);
    const auto norm_array(d["normalization"].GetArray());

    // A string of length L has at most L significant digits
    size_t digits10 = boost_float_thread_default_digits10();
    for(const auto &value : norm_array)
      digits10 = std::max<size_t>(digits10, value.GetStringLength());
    Scoped_Boost_Float_Precision precision(
      std::ceil(digits10 * std::log2(10.0)));

    std::vector<Boost_Float> normalization;
    normalization.reserve(norm_array.Size());
    for(const auto &value : norm_array)
      normalization.emplace_back(value.GetString());
    return normalization;
  }
}

std::optional<std::vector<Boost_Float>>
read_normalization(const fs::path &sdp_path, Timers &timers)
{
  Scoped_Timer timer(timers, "read_normalization");
  ASSERT(fs::exists(sdp_path), "SDP path does not exist: ", sdp_path);

  const std::string normalization_name("normalization.json");
  if(is_regular_file(sdp_path))
    {
      Archive_Reader reader(sdp_path);
      while(reader.next_entry())
        {
          if(normalization_name == archive_entry_pathname(reader.entry_ptr))
            {
              std::istream stream(&reader);
              return read_normalization_stream(stream);
            }
        }
    }
  else
    {
      auto normalization_path = sdp_path / normalization_name;
      if(exists(normalization_path))
        {
          std::ifstream normalization_stream(normalization_path);
          return read_normalization_stream(normalization_stream);
        }
    }
  return {};
}
