#pragma once

#include <El.hpp>

#include <archive.h>
#include <archive_entry.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <streambuf>
#include <memory>
#include <array>

struct Archive_Reader : public std::streambuf
{
  std::unique_ptr<archive, int (*)(archive *)> ptr;
  std::array<char, 8192> buffer;
  archive_entry *entry_ptr;
  bool entry_is_valid = false;
  Archive_Reader(const boost::filesystem::path &filename);
  Archive_Reader() = delete;
  ~Archive_Reader() = default;

  bool next_entry()
  {
    auto read_result = archive_read_next_header(ptr.get(), &entry_ptr);
    if(read_result == ARCHIVE_OK)
      entry_is_valid = true;
    else if(read_result == ARCHIVE_EOF)
      entry_is_valid = false;
    else
      El::RuntimeError("Archive_Reader: ", archive_error_string(ptr.get()));

    return entry_is_valid;
  }
  int underflow();
};
