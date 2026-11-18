#pragma once

#include "sdpb_util/Environment.hxx"
#include "sdpb_util/assert.hxx"

#include <El.hpp>

#include <cstddef>
#include <cstdlib>
#include <functional>
#include <vector>

// The copy of a small [MC,MR] block (e.g. a diagonal block A11 of
// Cholesky_Q, nb x nb) to [STAR,STAR] on every rank, with one copy per node
// instead of one per rank (--qCholeskyNodeDiag in the multimodular
// Cholesky_Q; --qSolveCache for the Q solve's diagonal blocks).
//
// El's copy (copy::AllGather) is one MPI_Allgather over all p ranks of
// the grid, so every rank receives all nb^2 serialized BigFloats and each
// node receives them (ranks per node) times. Here every rank serializes
// its local entries into a buffer in node-shared memory (an
// MPI_Win_allocate_shared window), one leader rank per node exchanges the
// nodes' parts with the other leaders (MPI_Allgatherv, in place, straight
// into the shared buffer), and every rank deserializes the whole block
// from the shared buffer. The bytes are El's serialization (BigFloat::
// Serialize / Deserialize), so the result is bit for bit El's copy.
//
// Layout of the shared buffer: the nodes' parts in node order (leaders in
// grid-rank order); inside a node's part, its ranks in VC (grid) order;
// each rank's local entries column by column, as El packs them.
// Two buffers alternate between calls: a rank can only write call k+2's
// entries after the node barrier of call k+1, i.e. after every rank of
// its node has finished reading call k.
//
// Nodes are Environment's nodes (MPI_COMM_TYPE_SHARED groups, or the groups
// of the test hook SDPB_FAKE_NODES, which never span machines), restricted to
// the grid's ranks.
//
// The same window also serves node-local exchanges (node_all_gather), used
// to share the factorization of a diagonal block among the node's ranks
// (--qCholeskyNodeFactor).
// Unit tests: the node index of this process for Node_Diag_Exchange, to try
// several splits of one host's ranks into nodes in one run (-1: Environment's
// node index).
inline int &node_diag_exchange_test_node()
{
  static int node = -1;
  return node;
}

class Node_Diag_Exchange
{
public:
  // max_entries: the largest matrix to be exchanged (height x width)
  Node_Diag_Exchange(const El::Grid &grid, const size_t max_entries)
      : grid(grid)
  {
    MPI_Comm vc_comm = grid.VCComm().comm;
    const int vc_rank = grid.VCRank(), size = grid.Size();
    // The grid's ranks on this process's node, as Environment defines nodes
    // (also with the test hook SDPB_FAKE_NODES)
    const int node = node_diag_exchange_test_node() >= 0
                       ? node_diag_exchange_test_node()
                       : Environment::process_node_index();
    ASSERT(node >= 0, "Node_Diag_Exchange needs an initialized Environment");
    MPI_Comm_split(vc_comm, node, vc_rank, &node_comm);
    MPI_Comm_rank(node_comm, &node_rank);
    MPI_Comm_size(node_comm, &node_size);
    MPI_Comm_split(vc_comm, node_rank == 0 ? 0 : MPI_UNDEFINED, vc_rank,
                   &leader_comm);
    int node_index = 0;
    if(node_rank == 0)
      {
        MPI_Comm_rank(leader_comm, &node_index);
        MPI_Comm_size(leader_comm, &num_nodes);
      }
    MPI_Bcast(&node_index, 1, MPI_INT, 0, node_comm);
    MPI_Bcast(&num_nodes, 1, MPI_INT, 0, node_comm);
    vc_node.resize(size);
    MPI_Allgather(&node_index, 1, MPI_INT, vc_node.data(), 1, MPI_INT,
                  vc_comm);
    // VC ranks in node order, VC order inside a node
    for(int node = 0; node < num_nodes; ++node)
      for(int vc = 0; vc < size; ++vc)
        if(vc_node[vc] == node)
          node_vcs.push_back(vc);
    ASSERT_EQUAL(node_vcs.size(), size_t(size));

    entry_bytes = El::BigFloat().SerializedSize();
    buffer_bytes = max_entries * entry_bytes;
    const MPI_Aint my_bytes = node_rank == 0 ? MPI_Aint(2 * buffer_bytes) : 0;
    std::byte *base = nullptr;
    MPI_Win_allocate_shared(my_bytes, 1, MPI_INFO_NULL, node_comm, &base,
                            &win);
    MPI_Aint bytes;
    int disp_unit;
    MPI_Win_shared_query(win, 0, &bytes, &disp_unit, &base);
    ASSERT_EQUAL(size_t(bytes), 2 * buffer_bytes);
    buffers[0] = base;
    buffers[1] = base + buffer_bytes;
    MPI_Win_lock_all(MPI_MODE_NOCHECK, win);
    offsets.resize(size);
    node_bytes.resize(num_nodes);
    node_displs.resize(num_nodes);
  }
  ~Node_Diag_Exchange()
  {
    MPI_Win_unlock_all(win);
    MPI_Win_free(&win);
    if(leader_comm != MPI_COMM_NULL)
      MPI_Comm_free(&leader_comm);
    MPI_Comm_free(&node_comm);
  }
  Node_Diag_Exchange(const Node_Diag_Exchange &) = delete;
  Node_Diag_Exchange &operator=(const Node_Diag_Exchange &) = delete;

  // B := A ([MC,MR] nb x nb, e.g. a view of the diagonal block) on every
  // rank of the grid, as El's B = A.
  void all_gather(const El::DistMatrix<El::BigFloat> &A,
                  El::DistMatrix<El::BigFloat, El::STAR, El::STAR> &B)
  {
    const El::Int height = A.Height(), width = A.Width();
    B.Resize(height, width);
    const int h = A.ColStride(), w = A.RowStride();
    const int col_align = A.ColAlign(), row_align = A.RowAlign();
    const std::byte *buf = exchange(
      [&](const int vc) {
        return size_t(local_length(height, vc % h, col_align, h)
                      * local_length(width, vc / h, row_align, w));
      },
      A.LockedMatrix());
    // unpack every rank's part
    auto &b = B.Matrix();
    for(int vc = 0; vc < grid.Size(); ++vc)
      {
        const int r = vc % h, c = vc / h;
        const int i0 = El::Shift(r, col_align, h),
                  j0 = El::Shift(c, row_align, w);
        const El::Int m = local_length(height, r, col_align, h),
                      n = local_length(width, c, row_align, w);
        const auto *p = reinterpret_cast<const El::byte *>(
          buf + offsets[vc] * entry_bytes);
        for(El::Int j = 0; j < n; ++j)
          for(El::Int i = 0; i < m; ++i)
            p = b.Ref(i0 + i * h, j0 + j * w).Deserialize(p);
      }
  }

  // Node-local all-gather (no exchange between nodes): node rank r's matrix
  // `mine` (count(r) entries, column by column) is put in the shared buffer
  // at node_offset(r) (in entries) for every rank of the node to read.
  // Used to split the diagonal block's factorization among a node's ranks.
  const std::byte *node_all_gather(const std::function<size_t(int)> &count,
                                   const El::Matrix<El::BigFloat> &mine)
  {
    node_offsets.resize(node_size + 1);
    node_offsets[0] = 0;
    for(int r = 0; r < node_size; ++r)
      node_offsets[r + 1] = node_offsets[r] + count(r);
    ASSERT(node_offsets[node_size] * entry_bytes <= buffer_bytes);
    ASSERT_EQUAL(size_t(mine.Height() * mine.Width()), count(node_rank));
    std::byte *buf = buffers[parity];
    parity ^= 1;
    auto *p
      = reinterpret_cast<El::byte *>(buf + node_offsets[node_rank] * entry_bytes);
    for(El::Int j = 0; j < mine.Width(); ++j)
      for(El::Int i = 0; i < mine.Height(); ++i)
        p = mine.CRef(i, j).Serialize(p);
    node_sync();
    return buf;
  }
  size_t node_offset(const int r) const { return node_offsets.at(r); }
  int ranks_on_node() const { return node_size; }
  int rank_on_node() const { return node_rank; }
  size_t entry_size() const { return entry_bytes; }

private:
  // Every rank's local matrix (count(vc) entries of VC rank vc, column by
  // column) into the shared buffer, then the nodes' parts to every node.
  // Returns the buffer; offsets[vc] is VC rank vc's part, in entries.
  const std::byte *exchange(const std::function<size_t(int)> &count,
                            const El::Matrix<El::BigFloat> &mine)
  {
    size_t pos = 0;
    int node = -1;
    for(const int vc : node_vcs)
      {
        if(vc_node[vc] != node)
          {
            if(node >= 0)
              node_bytes[node] = int((pos - node_displs[node]) * entry_bytes);
            node = vc_node[vc];
            node_displs[node] = int(pos);
          }
        offsets[vc] = pos;
        pos += count(vc);
      }
    node_bytes[node] = int((pos - node_displs[node]) * entry_bytes);
    for(auto &d : node_displs)
      d = int(d * entry_bytes);
    ASSERT(pos * entry_bytes <= buffer_bytes, DEBUG_STRING(pos),
           DEBUG_STRING(buffer_bytes));
    ASSERT_EQUAL(size_t(mine.Height() * mine.Width()),
                 count(grid.VCRank()));

    std::byte *buf = buffers[parity];
    parity ^= 1;
    // pack my entries, column by column
    auto *p = reinterpret_cast<El::byte *>(buf + offsets[grid.VCRank()]
                                           * entry_bytes);
    for(El::Int j = 0; j < mine.Width(); ++j)
      for(El::Int i = 0; i < mine.Height(); ++i)
        p = mine.CRef(i, j).Serialize(p);
    node_sync();
    if(leader_comm != MPI_COMM_NULL && num_nodes > 1)
      MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, buf,
                     node_bytes.data(), node_displs.data(), MPI_BYTE,
                     leader_comm);
    node_sync();
    return buf;
  }
  static El::Int
  local_length(const El::Int n, const int rank, const int align,
               const int stride)
  {
    return El::Length(n, El::Shift(rank, align, stride), stride);
  }
  void node_sync() const
  {
    MPI_Win_sync(win);
    MPI_Barrier(node_comm);
    MPI_Win_sync(win);
  }

  const El::Grid &grid;
  MPI_Comm node_comm = MPI_COMM_NULL, leader_comm = MPI_COMM_NULL;
  int node_rank = 0, node_size = 1, num_nodes = 1;
  std::vector<int> vc_node, node_vcs;
  std::vector<size_t> offsets, node_offsets;
  std::vector<int> node_bytes, node_displs;
  size_t entry_bytes = 0, buffer_bytes = 0;
  MPI_Win win = MPI_WIN_NULL;
  std::byte *buffers[2] = {nullptr, nullptr};
  int parity = 0;
};
