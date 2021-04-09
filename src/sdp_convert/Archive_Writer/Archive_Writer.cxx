#include "../Archive_Writer.hxx"

#include <stdexcept>

Archive_Writer::Archive_Writer(const boost::filesystem::path &filename)
    : ptr(archive_write_new(), archive_write_free)
{
  // Hard code zip with no compression.
  if(archive_write_set_format_zip(ptr.get()) != ARCHIVE_OK
     || archive_write_set_format_option(ptr.get(), "zip", "compression",
                                        "store")
          != ARCHIVE_OK
     || archive_write_open_filename(ptr.get(), filename.c_str()) != ARCHIVE_OK)
    {
      throw std::runtime_error(
        "Unable to set options for writing an archive.");
    }
}
