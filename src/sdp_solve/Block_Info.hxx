#pragma once

#include "sdpb_util/block_mapping/Abstract_Block_Info.hxx"

#include <filesystem>

class Block_Info final : public Abstract_Block_Info
{
public:
  std::vector<size_t> dimensions;
  std::vector<size_t> num_points;

  // TODO add Block_Info::total_size() == block_info.dimensions.size()
  // Rename block_info.block_indices to local_block_indices to make it clearer

  Block_Info() = delete;

  Block_Info(const Environment &env, const std::vector<size_t> &dimensions,
             const std::vector<size_t> &num_points,
             const std::vector<Block_Cost> &block_costs,
             const size_t &proc_granularity, const Verbosity &verbosity);

  static Block_Info
  create(const Environment &env, const std::filesystem::path &sdp_path,
         const std::filesystem::path &block_timings_path,
         const size_t &proc_granularity, const Verbosity &verbosity);
  static Block_Info
  create(const Environment &env, const std::filesystem::path &sdp_path,
         const El::Matrix<int32_t> &block_timings,
         const size_t &proc_granularity, const Verbosity &verbosity);
  static Block_Info
  create(const Environment &env, const std::vector<size_t> &dimensions,
         const std::vector<size_t> &num_points, const size_t &dual_dimension,
         const size_t &proc_granularity, const Verbosity &verbosity);

  [[nodiscard]] size_t get_schur_block_size(size_t index) const;
  [[nodiscard]] std::vector<size_t> schur_block_sizes() const;
  [[nodiscard]] size_t
  get_bilinear_pairing_block_size(size_t index, size_t parity) const;
  [[nodiscard]] std::vector<size_t> bilinear_pairing_block_sizes() const;
  [[nodiscard]] size_t
  get_psd_matrix_block_size(size_t index, size_t parity) const;
  [[nodiscard]] std::vector<size_t> psd_matrix_block_sizes() const;
  [[nodiscard]] size_t
  get_bilinear_bases_height(size_t index, size_t parity) const;
  [[nodiscard]] size_t
  get_bilinear_bases_width(size_t index, size_t /*parity*/) const;
};

void swap(Block_Info &a, Block_Info &b) noexcept;
