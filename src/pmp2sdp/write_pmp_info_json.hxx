#pragma once

#include "write_vector.hxx"
#include "pmp/PMP_Info.hxx"
#include "Reuse_SDP.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/json/Json_Writer.hxx"

#include <El.hpp>

#include <cstring>
#include <limits>
#include <ostream>
#include <type_traits>

inline void write_pvm_info_json(Json_Writer &writer, const PVM_Info &block,
                                const size_t num_blocks)
{
  block.validate(num_blocks);
  writer.StartObject();
  writer.Key("index");
  {
    writer.Int(block.block_index);
  }
  writer.Key("path");
  {
    writer.String(block.block_path.string().c_str());
  }
  writer.Key("dim");
  {
    writer.Uint64(block.dim);
  }
  writer.Key("prefactor");
  {
    writer.StartObject();
    writer.Key("constant");
    writer.BigFloat(block.prefactor.constant);
    writer.Key("base");
    writer.BigFloat(block.prefactor.base);
    writer.Key("poles");
    writer.BigFloatArray(block.prefactor.poles);
    writer.EndObject();
  }
  writer.Key("reducedPrefactor");
  {
    writer.StartObject();
    writer.Key("constant");
    writer.BigFloat(block.reduced_prefactor.constant);
    writer.Key("base");
    writer.BigFloat(block.reduced_prefactor.base);
    writer.Key("poles");
    writer.BigFloatArray(block.reduced_prefactor.poles);
    writer.EndObject();
  }
  writer.Key("samplePoints");
  writer.BigFloatArray(block.sample_points);
  writer.Key("sampleScalings");
  writer.BigFloatArray(block.sample_scalings);
  writer.Key("reducedSampleScalings");
  writer.BigFloatArray(block.reduced_sample_scalings);
  writer.EndObject();
}

inline void write_pmp_info_json(std::ostream &output_stream,
                                const std::vector<PVM_Info> &pmp_info)
{
  ASSERT_EQUAL(El::mpi::Rank(), 0);
  rapidjson::OStreamWrapper ows(output_stream);
  Json_Writer writer(ows);

  writer.StartArray();
  for(const auto &block : pmp_info)
    write_pvm_info_json(writer, block, pmp_info.size());

  writer.EndArray();
  ASSERT(writer.IsComplete());
}

// Same output as write_pmp_info_json(), for an SDP with blocks reused from
// an old SDP: reused entries are copied from the old pmp_info.json.
inline void write_pmp_info_json(std::ostream &output_stream,
                                const std::vector<PVM_Info> &pmp_info,
                                const Reused_SDP &reuse)
{
  ASSERT_EQUAL(El::mpi::Rank(), 0);
  ASSERT_EQUAL(reuse.new_to_old.size(), pmp_info.size());
  rapidjson::OStreamWrapper ows(output_stream);
  output_stream << '[';
  for(size_t index = 0; index < pmp_info.size(); ++index)
    {
      if(index != 0)
        output_stream << ',';
      const auto &old_index = reuse.new_to_old.at(index);
      if(old_index.has_value())
        {
          output_stream << reuse.renumbered_pmp_info_entry(old_index.value(),
                                                           index);
        }
      else
        {
          ows.Flush();
          Json_Writer writer(ows);
          write_pvm_info_json(writer, pmp_info.at(index), pmp_info.size());
          ASSERT(writer.IsComplete());
          ows.Flush();
        }
    }
  output_stream << ']';
}

// PVM_Info is collected on rank 0 in one MPI_Gatherv of a byte buffer per
// rank.
// Boost_Float values are packed exactly: their MPFR precision and an exact
// hex string. (Converting them to El::BigFloat and back rounds them, so
// pmp_info.json depended on which rank owned each block.)
class PVM_Info_Packer
{
public:
  std::vector<char> buffer;

  template <class T> void put(const T &value)
  {
    static_assert(std::is_trivially_copyable_v<T>);
    const auto *bytes = reinterpret_cast<const char *>(&value);
    buffer.insert(buffer.end(), bytes, bytes + sizeof(T));
  }
  void put_string(const std::string &str)
  {
    put<size_t>(str.size());
    buffer.insert(buffer.end(), str.begin(), str.end());
  }
  void put_boost_float(const Boost_Float &value)
  {
    const auto &x = value.backend().data();
    char *str = nullptr;
    ASSERT(mpfr_asprintf(&str, "%Ra", x) >= 0);
    const std::string hex(str);
    mpfr_free_str(str);
    put<size_t>(mpfr_get_prec(x));
    put_string(hex);
  }
  void put_big_floats(const std::vector<El::BigFloat> &vec)
  {
    put<size_t>(vec.size());
    for(const auto &value : vec)
      {
        const size_t offset = buffer.size();
        buffer.resize(offset + value.SerializedSize());
        value.Serialize(reinterpret_cast<El::byte *>(buffer.data() + offset));
      }
  }
  void put_damped_rational(const Damped_Rational &damped_rational)
  {
    put_boost_float(damped_rational.constant);
    put_boost_float(damped_rational.base);
    put<size_t>(damped_rational.poles.size());
    for(const auto &pole : damped_rational.poles)
      put_boost_float(pole);
  }
  void put_pvm_info(const PVM_Info &pvm_info)
  {
    put<int>(pvm_info.block_index);
    put_string(pvm_info.block_path.string());
    put<size_t>(pvm_info.dim);
    put_damped_rational(pvm_info.prefactor);
    put_damped_rational(pvm_info.reduced_prefactor);
    put_big_floats(pvm_info.sample_points);
    put_big_floats(pvm_info.sample_scalings);
    put_big_floats(pvm_info.reduced_sample_scalings);
  }
};

class PVM_Info_Unpacker
{
public:
  const char *pos;
  const char *const end;

  PVM_Info_Unpacker(const char *begin, const char *end)
      : pos(begin), end(end)
  {}

  template <class T> T get()
  {
    static_assert(std::is_trivially_copyable_v<T>);
    ASSERT(pos + sizeof(T) <= end, "PVM_Info buffer is too short");
    T value;
    std::memcpy(&value, pos, sizeof(T));
    pos += sizeof(T);
    return value;
  }
  std::string get_string()
  {
    const auto size = get<size_t>();
    ASSERT(pos + size <= end, "PVM_Info buffer is too short");
    std::string str(pos, size);
    pos += size;
    return str;
  }
  Boost_Float get_boost_float()
  {
    const auto prec = get<size_t>();
    const auto hex = get_string();
    Boost_Float result;
    auto &x = result.backend().data();
    mpfr_set_prec(x, prec);
    ASSERT(mpfr_set_str(x, hex.c_str(), 0, MPFR_RNDN) == 0, "Cannot parse ",
           hex);
    return result;
  }
  std::vector<El::BigFloat> get_big_floats()
  {
    std::vector<El::BigFloat> vec(get<size_t>());
    for(auto &value : vec)
      {
        ASSERT(pos + value.SerializedSize() <= end,
               "PVM_Info buffer is too short");
        pos = reinterpret_cast<const char *>(
          value.Deserialize(reinterpret_cast<const El::byte *>(pos)));
      }
    return vec;
  }
  Damped_Rational get_damped_rational()
  {
    Damped_Rational damped_rational;
    damped_rational.constant = get_boost_float();
    damped_rational.base = get_boost_float();
    damped_rational.poles.resize(get<size_t>());
    for(auto &pole : damped_rational.poles)
      pole = get_boost_float();
    return damped_rational;
  }
  PVM_Info get_pvm_info()
  {
    PVM_Info pvm_info;
    pvm_info.block_index = get<int>();
    pvm_info.block_path = get_string();
    pvm_info.dim = get<size_t>();
    pvm_info.prefactor = get_damped_rational();
    pvm_info.reduced_prefactor = get_damped_rational();
    pvm_info.sample_points = get_big_floats();
    pvm_info.sample_scalings = get_big_floats();
    pvm_info.reduced_sample_scalings = get_big_floats();
    return pvm_info;
  }
};

// Collect all blocks on rank=0.
// Blocks reused by pmp2sdp --reuseSdp are on no rank and stay empty.
inline std::vector<PVM_Info>
synchronize_pmp_info_blocks(const PMP_Info &pmp_info)
{
  const auto comm = El::mpi::COMM_WORLD;
  const int rank = comm.Rank();
  const int num_ranks = comm.Size();
  std::vector<PVM_Info> all_pvm_infos(pmp_info.num_blocks);

  PVM_Info_Packer packer;
  if(rank != 0)
    {
      packer.put<size_t>(pmp_info.blocks.size());
      for(const auto &block : pmp_info.blocks)
        packer.put_pvm_info(block);
    }
  // MPI counts are int
  ASSERT(packer.buffer.size()
           <= static_cast<size_t>(std::numeric_limits<int>::max()),
         DEBUG_STRING(packer.buffer.size()));
  const int send_size = packer.buffer.size();

  std::vector<int> recv_sizes(rank == 0 ? num_ranks : 0);
  MPI_Gather(&send_size, 1, MPI_INT, recv_sizes.data(), 1, MPI_INT, 0,
             comm.comm);

  std::vector<int> displacements;
  std::vector<char> recv_buffer;
  if(rank == 0)
    {
      size_t total_size = 0;
      for(const int size : recv_sizes)
        {
          displacements.push_back(total_size);
          total_size += size;
          ASSERT(total_size
                   <= static_cast<size_t>(std::numeric_limits<int>::max()),
                 DEBUG_STRING(total_size));
        }
      recv_buffer.resize(total_size);
    }
  MPI_Gatherv(packer.buffer.data(), send_size, MPI_CHAR, recv_buffer.data(),
              recv_sizes.data(), displacements.data(), MPI_CHAR, 0,
              comm.comm);

  const auto set_block = [&](PVM_Info &&block) {
    const int block_index = block.block_index;
    ASSERT(block_index >= 0 && block_index < pmp_info.num_blocks,
           DEBUG_STRING(block_index), DEBUG_STRING(pmp_info.num_blocks));
    ASSERT(all_pvm_infos.at(block_index).block_index == -1,
           "block_", block_index, " is on more than one rank");
    all_pvm_infos.at(block_index) = std::move(block);
  };
  // Every rank keeps its own blocks, as before
  for(const auto &block : pmp_info.blocks)
    set_block(PVM_Info(block));
  if(rank != 0)
    return all_pvm_infos;
  for(int from = 1; from < num_ranks; ++from)
    {
      const char *begin = recv_buffer.data() + displacements.at(from);
      PVM_Info_Unpacker unpacker(begin, begin + recv_sizes.at(from));
      const auto num_blocks = unpacker.get<size_t>();
      for(size_t i = 0; i < num_blocks; ++i)
        set_block(unpacker.get_pvm_info());
      ASSERT(unpacker.pos == unpacker.end, DEBUG_STRING(from));
    }
  return all_pvm_infos;
}
