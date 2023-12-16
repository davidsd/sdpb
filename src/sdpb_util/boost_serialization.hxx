#pragma once

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/version.hpp>

#include <El.hpp>

namespace boost::serialization
{
  // El::BigFloat

  template <class Archive>
  void serialize(Archive &ar, El::BigFloat &f,
                 const boost::serialization::version_type &version)
  {
    // SDP blocks may contain many zeros.
    // In that case binary format can be even less compact that json,
    // see e.g. https://github.com/davidsd/sdpb/issues/148
    // To optimize storing zeros, we use boolean flag is_zero.
    // if (is_zero == true), we don't serialize all BigFloat bytes.
    // Optimization introduced since version = 1

    const El::BigFloat zero(0);
    bool is_zero = false;
    if(version > 0)
      {
        if(Archive::is_saving::value)
          is_zero = (f == zero);
        ar & is_zero;
      }

    if(is_zero)
      {
        if(Archive::is_loading::value)
          f = zero;
        return;
      }

    auto size = f.SerializedSize();
    std::vector<El::byte> vec(size);
    El::byte *buffer = vec.data();
    // array serialization is more compact than vector;
    // for precision=512 the difference is ~9%.
    auto array = make_array(buffer, size);

    if(Archive::is_saving::value)
      f.Serialize(buffer);

    ar & array;

    if(Archive::is_loading::value)
      f.Deserialize(buffer);
  }

  // El::Matrix

  template <class Archive, class Ring>
  void save(Archive &ar, El::Matrix<Ring> const &matrix,
            const boost::serialization::version_type &)
  {
    ar & matrix.Height();
    ar & matrix.Width();
    ar & matrix.LDim();
    auto size = matrix.LDim() * matrix.Width();
    auto data = boost::serialization::make_array(matrix.LockedBuffer(), size);
    ar & data;
  }
  template <class Archive, class Ring>
  void load(Archive &ar, El::Matrix<Ring> &matrix,
            const boost::serialization::version_type &)
  {
    El::Int height, width, leadingDimension;
    ar & height;
    ar & width;
    ar & leadingDimension;
    matrix.Resize(height, width, leadingDimension);
    auto size = leadingDimension * width;
    auto data = boost::serialization::make_array(matrix.Buffer(), size);
    ar & data;
  }
}

BOOST_CLASS_VERSION(El::BigFloat, 1)

BOOST_SERIALIZATION_SPLIT_FREE(El::Matrix<El::BigFloat>)

// https://www.boost.org/doc/libs/1_82_0/libs/serialization/doc/special.html#objecttracking
// We are just writing arrays, don't need the object tracking mechanism,
// which may cause memory issues (according to some StackOverflow questions).
BOOST_CLASS_TRACKING(El::BigFloat, boost::serialization::track_never)
BOOST_CLASS_TRACKING(El::Matrix<El::BigFloat>,
                     boost::serialization::track_never)
