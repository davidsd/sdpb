#include "Archive_Entry.hxx"

Archive_Entry::Archive_Entry(const boost::filesystem::path &filename,
                             const int64_t &num_bytes)
    : entry_ptr(archive_entry_new(), archive_entry_free)
{
  archive_entry_set_pathname(entry_ptr.get(), filename.filename().c_str());
  archive_entry_set_size(entry_ptr.get(), num_bytes);
  archive_entry_set_filetype(entry_ptr.get(), AE_IFREG);
  archive_entry_set_perm(entry_ptr.get(), 0644);
}
