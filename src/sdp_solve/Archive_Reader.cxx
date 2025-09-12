#include "Archive_Reader.hxx"

#include <stdexcept>
#include <string>

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

int Archive_Reader::underflow()
{
  if(gptr() == egptr())
    {
      ssize_t num_bytes_read(
        archive_read_data(ptr.get(), buffer.data(), buffer.size()));
      if(num_bytes_read < 0)
        {
          RUNTIME_ERROR("Error reading archive.");
        }
      if(num_bytes_read > 0)
        {
          setg(buffer.data(), buffer.data(), buffer.data() + num_bytes_read);
        }
    }
  return (gptr() == egptr()) ? std::char_traits<char>::eof()
                             : std::char_traits<char>::to_int_type(*gptr());
}
