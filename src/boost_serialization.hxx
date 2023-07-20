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
                 const boost::serialization::version_type &)
  {
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
    auto size = leadingDimension * width;
    Ring *buffer = new Ring[size];
    auto data = boost::serialization::make_array(buffer, size);
    ar & data;
    matrix.Control(height, width, buffer, leadingDimension);
  }
}

BOOST_SERIALIZATION_SPLIT_FREE(El::Matrix<El::BigFloat>)
