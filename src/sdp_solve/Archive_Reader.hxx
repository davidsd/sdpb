#pragma once

#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <archive.h>
#include <archive_entry.h>

#include <filesystem>
#include <streambuf>
#include <memory>
#include <array>

struct Archive_Reader : public std::streambuf
{
  std::unique_ptr<archive, int (*)(archive *)> ptr;
  std::array<char, 8192> buffer;
  archive_entry *entry_ptr;
  bool entry_is_valid = false;
  explicit Archive_Reader(const std::filesystem::path &filename);
  Archive_Reader() = delete;
  ~Archive_Reader() override = default;

  bool next_entry()
  {
    const auto read_result = archive_read_next_header(ptr.get(), &entry_ptr);
    if(read_result == ARCHIVE_OK)
      entry_is_valid = true;
    else if(read_result == ARCHIVE_EOF)
      entry_is_valid = false;
    else
      {
        const char *const pathname = archive_entry_pathname(entry_ptr);
        const std::string name = pathname ? pathname : "EOF";
        RUNTIME_ERROR("Archive_Reader: ", name, ": ",
                         archive_error_string(ptr.get()));
      }

    return entry_is_valid;
  }
  int underflow() override;
};
