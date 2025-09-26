#include "Archive_Reader.hxx"

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

bool Archive_Reader::next_entry()
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
