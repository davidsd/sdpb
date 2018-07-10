#pragma once

#include <libxml++/libxml++.h>

#include <cstring>
#include <string>

// The default Glib::ustring==std::string converts everything to UCS4
// and is terribly slow.  This only checks for equality and only works
// if std_string is pure ASCII.
inline bool glib_equals_string(const Glib::ustring &glib_string,
                               const std::string &std_string)
{
  return glib_string.bytes() == std_string.size()
         && std::strcmp(glib_string.data(), std_string.data()) == 0;
}
