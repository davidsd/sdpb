#pragma once

#include "write_vector.hxx"
#include "pmp/Damped_Rational.hxx"
#include "pmp/Polynomial_Matrix_Program.hxx"
#include "pmp/PVM_Info.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Boost_Float.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"

#include <El.hpp>

#include <iostream>
#include <filesystem>
#include <optional>
#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/writer.h>

namespace fs = std::filesystem;

inline void write_pmp_info_json(std::ostream &output_stream,
                                const std::vector<PVM_Info> &pmp_info)
{
  ASSERT_EQUAL(El::mpi::Rank(), 0);
  const size_t num_blocks = pmp_info.size();
  rapidjson::OStreamWrapper ows(output_stream);
  rapidjson::Writer writer(ows);

  // reusable stream for big floats
  std::stringstream ss;
  set_stream_precision(ss);
  auto add_BigFloat = [&writer, &ss](const El::BigFloat &value) {
    ss.str({});
    ss << value;
    writer.String(ss.str().c_str());
  };
  auto add_Boost_Float = [&writer, &ss](const Boost_Float &value) {
    ss.str({});
    ss << value;
    writer.String(ss.str().c_str());
  };

  auto add_bigFloat_array = [&](const std::vector<El::BigFloat> &value) {
    writer.StartArray();
    for(auto &p : value)
      add_BigFloat(p);
    writer.EndArray();
  };
  auto add_Boost_Float_array = [&](const std::vector<Boost_Float> &value) {
    writer.StartArray();
    for(auto &p : value)
      add_Boost_Float(p);
    writer.EndArray();
  };

  writer.StartArray();
  for(const auto &block : pmp_info)
    {
      writer.StartObject();
      writer.Key("block_path");
      {
        writer.String(block.block_path.string().c_str());
      }
      writer.Key("prefactor");
      {
        writer.StartObject();
        writer.Key("constant");
        add_Boost_Float(block.prefactor.constant);
        writer.Key("base");
        add_Boost_Float(block.prefactor.base);
        writer.Key("poles");
        add_Boost_Float_array(block.prefactor.poles);
        writer.EndObject();
      }
      writer.Key("reducedPrefactor");
      {
        writer.StartObject();
        writer.Key("constant");
        add_Boost_Float(block.reduced_prefactor.constant);
        writer.Key("base");
        add_Boost_Float(block.reduced_prefactor.base);
        writer.Key("poles");
        add_Boost_Float_array(block.reduced_prefactor.poles);
        writer.EndObject();
      }
      writer.Key("samplePoints");
      add_bigFloat_array(block.sample_points);
      writer.Key("sampleScalings");
      add_bigFloat_array(block.sample_scalings);
      writer.Key("reducedSampleScalings");
      add_bigFloat_array(block.reduced_sample_scalings);
      writer.EndObject();
    }

  writer.EndArray();
  ASSERT(writer.IsComplete());
}

// Send from rank='from' to rank=0
inline void synchronize_pvm_info(PVM_Info &pvm_info, const int from)
{
  const int to = 0;

  const auto comm = El::mpi::COMM_WORLD;
  const int rank = comm.Rank();

  if(rank != to && rank != from)
    return;
  if(to == from)
    return;

  // block_path
  {
    auto block_path = pvm_info.block_path.string();
    if(rank == from)
      {
        const size_t path_size = block_path.size();
        El::mpi::Send<size_t>(path_size, to, comm);
        MPI_Send(block_path.data(), path_size, MPI_CHAR, to, 0, comm.comm);
      }
    if(rank == to)
      {
        const auto path_size = El::mpi::Recv<size_t>(from, comm);
        block_path.resize(path_size);
        MPI_Recv(block_path.data(), path_size, MPI_CHAR, from, 0, comm.comm,
                 MPI_STATUS_IGNORE);
        pvm_info.block_path = block_path;
      }
  }

  // prefactor, reduced_prefactor
  for(auto *damped_rational_ptr :
      {&pvm_info.prefactor, &pvm_info.reduced_prefactor})
    {
      auto &damped_rational = *damped_rational_ptr;
      if(rank == from)
        {
          El::mpi::Send(to_BigFloat(damped_rational.constant), to, comm);
          El::mpi::Send(to_BigFloat(damped_rational.base), to, comm);
          const auto poles = to_BigFloat_Vector(damped_rational.poles);
          El::mpi::Send<size_t>(poles.size(), to, comm);
          El::mpi::Send(poles.data(), poles.size(), to, comm);
        }
      if(rank == to)
        {
          damped_rational.constant
            = to_Boost_Float(El::mpi::Recv<El::BigFloat>(from, comm));
          damped_rational.base
            = to_Boost_Float(El::mpi::Recv<El::BigFloat>(from, comm));

          const size_t num_poles = El::mpi::Recv<size_t>(from, comm);
          std::vector<El::BigFloat> poles(num_poles);
          El::mpi::Recv(poles.data(), poles.size(), from, comm);
          damped_rational.poles = to_Boost_Float_Vector(poles);
        }
    }
  for(auto vec_ptr : {&pvm_info.sample_points, &pvm_info.sample_scalings,
                      &pvm_info.reduced_sample_scalings})
    {
      auto &vec = *vec_ptr;
      if(rank == from)
        {
          El::mpi::Send<size_t>(vec.size(), to, comm);
          El::mpi::Send(vec.data(), vec.size(), to, comm);
        }
      if(rank == to)
        {
          const size_t num_points = El::mpi::Recv<size_t>(from, comm);
          vec.resize(num_points);
          El::mpi::Recv(vec.data(), vec.size(), from, comm);
        }
    }
}

std::vector<PVM_Info> inline synchronize_pmp_info(
  const Polynomial_Matrix_Program &pmp)
{
  std::vector<PVM_Info> pmp_info(pmp.num_matrices);

  // Fill pmp_info with all local blocks
  for(size_t local_index = 0;
      local_index < pmp.matrix_index_local_to_global.size(); ++local_index)
    {
      const size_t block_index
        = pmp.matrix_index_local_to_global.at(local_index);
      const auto &pvm = pmp.matrices.at(local_index);

      auto &pvm_info = pmp_info.at(block_index);
      pvm_info.block_path = pmp.block_paths.at(local_index);
      pvm_info.prefactor = pvm.prefactor;
      pvm_info.reduced_prefactor = pvm.reduced_prefactor;
      pvm_info.sample_points = pvm.sample_points;
      pvm_info.sample_scalings = pvm.sample_scalings;
      pvm_info.reduced_sample_scalings = pvm.reduced_sample_scalings;
    }

  // Synchronize owning ranks for each block
  const int rank = El::mpi::Rank();
  std::vector<int> block_ranks(pmp.num_matrices, -1);
  for(auto &block_index : pmp.matrix_index_local_to_global)
    {
      block_ranks.at(block_index) = rank;
    }
  El::mpi::AllReduce(block_ranks.data(), block_ranks.size(), El::mpi::MAX,
                     El::mpi::COMM_WORLD);

  // Send all blocks to rank=0
  for(size_t block_index = 0; block_index < pmp.num_matrices; ++block_index)
    {
      const int from = block_ranks.at(block_index);
      synchronize_pvm_info(pmp_info.at(block_index), from);
    }
  return pmp_info;
}