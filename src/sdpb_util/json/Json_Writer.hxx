#pragma once

#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/writer.h>

#include <sstream>
#include <utility>

// Helper class to simplify writing JSON
// TBaseWriter is rapidjson::Writer<...> or rapidjson::PrettyWriter<...>
// TODO move to sdpb_util, reuse for write_pmp_info and save_c_minus_By
template <class TBaseWriter> class Json_BigFloat_Writer : public TBaseWriter
{
public:
  template <class... TArgs>
  explicit Json_BigFloat_Writer(TArgs &&...args)
      : TBaseWriter(std::forward<TArgs>(args)...)
  { set_stream_precision(ss); }

  auto BigFloat(const El::BigFloat &value) { return BigFloat_impl(value); }
  auto BigFloat(const Boost_Float &value) { return BigFloat_impl(value); }

  template <class TFloat> auto BigFloatArray(const std::vector<TFloat> &arr)
  {
    this->StartArray();
    for(const auto &value : arr)
      BigFloat(value);
    return this->EndArray();
  }

private:
  // Reusable stream for writing BigFloats
  std::stringstream ss;

  template <class T> auto BigFloat_impl(const T &value)
  {
    ss.str({});
    ss << value;
    return this->String(ss.str().c_str());
  }
};

using Json_Writer
  = Json_BigFloat_Writer<rapidjson::Writer<rapidjson::OStreamWrapper>>;
using Json_PrettyWriter
  = Json_BigFloat_Writer<rapidjson::PrettyWriter<rapidjson::OStreamWrapper>>;
