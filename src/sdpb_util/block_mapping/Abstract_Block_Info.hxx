#pragma once

#include "sdpb_util/Environment.hxx"
#include "sdpb_util/Verbosity.hxx"
#include "sdpb_util/block_mapping/Block_Cost.hxx"
#include "sdpb_util/block_mapping/Block_Mapping.hxx"
#include "sdpb_util/block_mapping/MPI_Comm_Wrapper.hxx"
#include "sdpb_util/block_mapping/MPI_Group_Wrapper.hxx"

#include <El.hpp>
#include <filesystem>

struct Abstract_Block_Info
{
  std::vector<size_t> block_indices;
  MPI_Group_Wrapper mpi_group;
  MPI_Comm_Wrapper mpi_comm;
  Block_Mapping block_mapping;

  // NB: Block_Info owns mpi_group and mpi_comm,
  // and clears them in destructor.
  // But it does not own node_comm, which comes from oustide (Environment).
  // Thus, we don't use MPI_Comm_Wrapper here.
  // TODO: refactor code to make semantics clearer.
  El::mpi::Comm node_comm;
  size_t node_index;

  // Index of current MPI group on a node
  [[nodiscard]] size_t node_group_index() const;
  [[nodiscard]] size_t num_nodes() const;

protected:
  Abstract_Block_Info() = delete;

  Abstract_Block_Info(
    const Environment &env, const std::vector<Block_Cost> &block_costs,
    const size_t &proc_granularity,
    const std::function<std::ostream &(std::ostream &, size_t)> &print_block,
    const Verbosity &verbosity);

  virtual ~Abstract_Block_Info() = default;
};

void swap(Abstract_Block_Info &a, Abstract_Block_Info &b) noexcept;
