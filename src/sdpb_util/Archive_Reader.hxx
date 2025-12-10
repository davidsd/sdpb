#pragma once

#include "assert.hxx"

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

  bool next_entry();

  int underflow() override;
};
