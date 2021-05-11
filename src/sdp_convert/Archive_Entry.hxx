#pragma once

#include <archive_entry.h>

#include <boost/filesystem.hpp>

#include <memory>
#include <array>
#include <stdexcept>

struct Archive_Entry
{
  std::unique_ptr<archive_entry, void (*)(archive_entry *)> entry_ptr;
  Archive_Entry(const boost::filesystem::path &filename,
                const int64_t &num_bytes);
};
