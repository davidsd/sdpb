#include "../Archive_Reader.hxx"

#include <stdexcept>

Archive_Reader::Archive_Reader(const std::filesystem::path &filename)
    : ptr(archive_read_new(), archive_read_free)
{
  if(archive_read_support_filter_all(ptr.get()) != ARCHIVE_OK
     || archive_read_support_format_all(ptr.get()) != ARCHIVE_OK
     || archive_read_open_filename(ptr.get(), filename.c_str(), 10240)
          != ARCHIVE_OK)
    {
      RUNTIME_ERROR("Error when opening ", filename);
    }
  setg(buffer.data(), buffer.data(), buffer.data());
}
