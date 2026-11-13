#pragma once

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/array.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include <El.hpp>

#include <algorithm>
#include <cstdlib>
#include <cstring>

// pmp2sdp --encodingPoints computes at a higher precision than the SDP is
// written at.  El::BigFloat::Serialize() writes gmp::num_limbs limbs, set by
// El's global precision (which can be set only once), and SDPB reads as many
// limbs as its own precision gives.  If output_num_limbs() > 0, BigFloats are
// saved in El's layout (int _mp_size, mp_exp_t _mp_exp, limbs) with that many
// limbs instead: the most significant ones, i.e. the value truncated.
// Loading is unchanged.  0 (default): El's own Serialize().
namespace sdpb_serialization
{
  inline int &output_num_limbs()
  {
    static int num_limbs = 0;
    return num_limbs;
  }
  // El's gmp::num_limbs for a precision (El's gmp.cpp, SetPrecision)
  inline int num_limbs_for_precision(const mp_bitcnt_t precision)
  {
    return (std::max(precision, static_cast<mp_bitcnt_t>(53))
            + 2 * GMP_NUMB_BITS - 1)
             / GMP_NUMB_BITS
           + 1;
  }
}

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

    const int output_num_limbs = sdpb_serialization::output_num_limbs();
    if(Archive::is_saving::value && output_num_limbs > 0)
      {
        const auto &mpf = f.gmp_float.get_mpf_t()[0];
        const int abs_size = std::abs(mpf._mp_size);
        const int n = std::min(abs_size, output_num_limbs);
        const int sign_size = mpf._mp_size < 0 ? -n : n;
        std::vector<El::byte> vec(sizeof(int) + sizeof(mp_exp_t)
                                  + output_num_limbs * sizeof(mp_limb_t));
        El::byte *buffer = vec.data();
        std::memcpy(buffer, &sign_size, sizeof(int));
        std::memcpy(buffer + sizeof(int), &mpf._mp_exp, sizeof(mp_exp_t));
        std::memcpy(buffer + sizeof(int) + sizeof(mp_exp_t),
                    mpf._mp_d + (abs_size - n), n * sizeof(mp_limb_t));
        auto array = make_array(buffer, vec.size());
        ar & array;
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
