#pragma once

#include <El.hpp>

#include <rapidjson/rapidjson.h>
#include <rapidjson/encodings.h>

#include <cstdint>
#include <rapidjson/reader.h>

// This abstract class duplicates rapidjson::BaseReaderHandler interface.
// We cannot store, e.g., vector of pointers to rapidjson::BaseReaderHandler<Encoding, Derived>
// because we have to specify the template parameter Derived (which is different for different classes)
// But we can store pointers to Abstract_Json_Reader_Handler instead.
struct Abstract_Json_Reader_Handler
    : rapidjson::BaseReaderHandler<rapidjson::UTF8<>,
                                   Abstract_Json_Reader_Handler>
{
  using SizeType = rapidjson::SizeType;
  using Ch = rapidjson::UTF8<>::Ch;

  // BaseReaderHandler interface
  // NB: do not ovverride these functions,
  // override virtual json_XXX() instead, see below!

  bool Default() { return json_default(); }
  bool Null() { return json_null(); }
  bool Bool(bool b) { return json_bool(b); }
  bool Int(int i) { return json_int(i); }
  bool Uint(unsigned i) { return json_uint(i); }
  bool Int64(int64_t i) { return json_int64(i); }
  bool Uint64(uint64_t i) { return json_uint64(i); }
  bool Double(double d) { return json_double(d); }
  bool RawNumber(const Ch *str, SizeType length, bool copy)
  {
    return json_raw_number(str, length, copy);
  }
  bool String(const Ch *str, SizeType length, bool copy)
  {
    return json_string(str, length, copy);
  }
  bool StartObject() { return json_start_object(); }
  bool Key(const Ch *str, SizeType length, bool copy)
  {
    return json_key(str, length, copy);
  }
  bool EndObject(SizeType memberCount) { return json_end_object(memberCount); }
  bool StartArray() { return json_start_array(); }
  bool EndArray(SizeType elementCount) { return json_end_array(elementCount); }

  // Implementing BaseReaderHandler interface via virtual functions

#define VIRTUAL_NOT_IMPLEMENTED(func)                                         \
  virtual func [[noreturn]]                                                   \
  {                                                                           \
    RUNTIME_ERROR("Not implemented: function '", #func, "' in class: '",      \
                  typeid(*this).name(), "'");                                 \
  }

  VIRTUAL_NOT_IMPLEMENTED(bool json_default())
  VIRTUAL_NOT_IMPLEMENTED(bool json_null())
  VIRTUAL_NOT_IMPLEMENTED(bool json_bool(bool b))
  VIRTUAL_NOT_IMPLEMENTED(bool json_int(int i))
  VIRTUAL_NOT_IMPLEMENTED(bool json_uint(unsigned i))
  VIRTUAL_NOT_IMPLEMENTED(bool json_int64(int64_t i))
  VIRTUAL_NOT_IMPLEMENTED(bool json_uint64(uint64_t i))
  VIRTUAL_NOT_IMPLEMENTED(bool json_double(double d))
  VIRTUAL_NOT_IMPLEMENTED(bool json_raw_number(const Ch *str, SizeType length,
                                               bool copy))
  VIRTUAL_NOT_IMPLEMENTED(bool json_string(const Ch *str, SizeType length,
                                           bool copy))
  VIRTUAL_NOT_IMPLEMENTED(bool json_start_object())
  VIRTUAL_NOT_IMPLEMENTED(bool json_key(const Ch *str, SizeType length,
                                        bool copy))
  VIRTUAL_NOT_IMPLEMENTED(bool json_end_object(SizeType memberCount))
  VIRTUAL_NOT_IMPLEMENTED(bool json_start_array())
  VIRTUAL_NOT_IMPLEMENTED(bool json_end_array(SizeType elementCount))

protected:
  ~Abstract_Json_Reader_Handler() = default;
};
