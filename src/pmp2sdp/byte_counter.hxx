// (C) Copyright 2008 CodeRage, LLC (turkanis at coderage dot com)
// (C) Copyright 2005-2007 Jonathan Turkanis
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt.)

// Copied from boost/include/boost/iostreams/counter.hpp
// Modified to
// 1) Use std::streamsize everywhere instead of int
// 2) Only count bytes, not lines
// 3) Only implements write, not read
// 4) Only chars, not wchars

#pragma once

#include <algorithm>
#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/char_traits.hpp>
#include <boost/iostreams/operations.hpp>
#include <boost/iostreams/pipeline.hpp>

class byte_counter
{
public:
  typedef char char_type;
  struct category : boost::iostreams::dual_use,
                    boost::iostreams::filter_tag,
                    boost::iostreams::multichar_tag,
                    boost::iostreams::optimally_buffered_tag
  {};
  explicit byte_counter(int first_byte = 0) : num_bytes(first_byte) {}
  std::streamsize optimal_buffer_size() const { return 0; }

  template <typename Sink>
  std::streamsize write(Sink &snk, const char *s, std::streamsize n)
  {
    std::streamsize result = boost::iostreams::write(snk, s, n);
    num_bytes += result;
    return result;
  }
  template <typename Source>
  std::streamsize read(Source &, char_type *, std::streamsize)
  {
    throw std::runtime_error(
      "INTERNAL_ERROR: byte_counter::read() not implemented");
  }

  std::streamsize num_bytes;
};
BOOST_IOSTREAMS_PIPABLE(byte_counter, 0)
