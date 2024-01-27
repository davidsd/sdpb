#include "../Archive_Reader.hxx"

#include <stdexcept>
#include <string>

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
