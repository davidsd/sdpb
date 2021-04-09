#include "../Archive_Writer.hxx"

#include <stdexcept>

#include <iostream>

void Archive_Writer::write_entry(const Archive_Entry &entry,
                                 std::istream &stream)
{
  if(archive_write_header(ptr.get(), entry.entry_ptr.get()) != ARCHIVE_OK)
    {
      throw std::runtime_error("Error when writing archive header");
    }
  std::array<char, 8192> buffer;
  stream.read(buffer.data(), buffer.size());
  while(stream.good())
    {
      if(archive_write_data(ptr.get(), buffer.data(), buffer.size()) < 0)
        {
          throw std::runtime_error("Error when writing archive data");
        }
      stream.read(buffer.data(), buffer.size());
    }
  if(archive_write_data(ptr.get(), buffer.data(), stream.gcount()) < 0)
    {
      throw std::runtime_error("Error when writing archive final data");
    }
}
