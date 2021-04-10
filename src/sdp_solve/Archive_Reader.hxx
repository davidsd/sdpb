#pragma once

#include <archive.h>
#include <archive_entry.h>

#include <boost/filesystem.hpp>

#include <streambuf>
#include <memory>
#include <array>

struct Archive_Reader : public std::streambuf
{
  std::unique_ptr<archive, int (*)(archive *)> ptr;
  std::array<char, 8192> buffer;
  archive_entry *entry_ptr;
  bool entry_is_valid  = false;
  Archive_Reader(const boost::filesystem::path &filename);
  Archive_Reader() = delete;
  ~Archive_Reader() = default;

  bool next_entry()
  {
    entry_is_valid=(archive_read_next_header(ptr.get(), &entry_ptr)
                    == ARCHIVE_OK);
    return entry_is_valid;
  }
  int underflow();
};
