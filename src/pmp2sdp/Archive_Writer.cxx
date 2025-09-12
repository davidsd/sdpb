#include "Archive_Writer.hxx"

#include "sdpb_util/assert.hxx"

#include <stdexcept>

Archive_Writer::Archive_Writer(const std::filesystem::path &filename)
    : ptr(archive_write_new(), archive_write_free)
{
  // Hard code zip with no compression.
  if(archive_write_set_format_zip(ptr.get()) != ARCHIVE_OK
     || archive_write_set_format_option(ptr.get(), "zip", "compression",
                                        "store")
          != ARCHIVE_OK
     || archive_write_open_filename(ptr.get(), filename.c_str()) != ARCHIVE_OK)
    {
      RUNTIME_ERROR("Unable to set options for writing an archive.");
    }
}

void Archive_Writer::write_entry(const Archive_Entry &entry,
                                 std::istream &stream)
{
  if(archive_write_header(ptr.get(), entry.entry_ptr.get()) != ARCHIVE_OK)
    {
      RUNTIME_ERROR("Error when writing archive header");
    }
  std::array<char, 8192> buffer;
  stream.read(buffer.data(), buffer.size());
  while(stream.good())
    {
      if(archive_write_data(ptr.get(), buffer.data(), buffer.size()) < 0)
        {
          RUNTIME_ERROR("Error when writing archive data");
        }
      stream.read(buffer.data(), buffer.size());
    }
  if(archive_write_data(ptr.get(), buffer.data(), stream.gcount()) < 0)
    {
      RUNTIME_ERROR("Error when writing archive final data");
    }
}
