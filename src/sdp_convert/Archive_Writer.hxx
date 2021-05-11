#pragma once

#include "Archive_Entry.hxx"
#include <archive.h>

#include <boost/filesystem.hpp>

#include <memory>

struct Archive_Writer
{
  std::unique_ptr<archive, int (*)(archive *)> ptr;
  Archive_Writer(const boost::filesystem::path &filename);
  ~Archive_Writer() { archive_write_close(ptr.get()); }
  Archive_Writer() = delete;
  
  void write_entry(const Archive_Entry &entry, std::istream &stream);
};
