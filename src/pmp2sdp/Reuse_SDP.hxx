#pragma once

#include "Block_File_Format.hxx"

#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <vector>

// Blocks of an existing SDP directory, to be reused by pmp2sdp --reuseSdp.
//
// An input file whose path appears in old_sdp/pmp_info.json is assumed
// unchanged since old_sdp was written: it is not parsed, and its blocks
// (block_info_*.json, block_data_*) are hard-linked (or copied) from old_sdp,
// renumbered to their new block indices. Its pmp_info.json entries are copied
// as raw text, with only "index" changed.
// All other input files are parsed and converted as usual.
//
// Requirements, checked at run time:
// - old_sdp is a directory written by pmp2sdp (not a zip),
//   with the same --precision, --maxNumPoles and --outputFormat;
// - the normalization vector read from the parsed files equals the old one
//   (normalization.json byte-equal, or absent in both),
//   because every block depends on it;
// - the objective is found in the parsed files, and has the old length.
// Not checked: that a reused file's content did not change (files are
// matched by path only), and that a reused file does not carry an
// objective different from the parsed one.
// Value of a pmp2sdp option in a command line as control.json records it
// (space-separated), or nullopt if absent. Reads "--name v", "--name=v",
// "-s v" and "-sv"; refuses an abbreviation of name ("--maxNum 2").
std::optional<std::string>
find_pmp2sdp_option(const std::vector<std::string> &tokens,
                    const std::string &name, const std::string &short_name);

struct Reused_SDP
{
  std::filesystem::path old_sdp;
  // Raw text of each old pmp_info.json entry, by old block index
  std::vector<std::string> raw_pmp_info_entries;
  // Old block indices of each input path, in increasing order
  std::map<std::filesystem::path, std::vector<size_t>> old_blocks_by_path;
  Block_File_Format format;

  // For each new block index: the old block index, if the block is reused.
  // Filled by set_new_block_mapping().
  std::vector<std::optional<size_t>> new_to_old;

  // Read on rank 0 and broadcast to all ranks.
  Reused_SDP(const std::filesystem::path &old_sdp, int precision,
             int64_t max_num_poles, Block_File_Format format);

  // Number of matrices of a reused input file, or nullopt if it must be parsed
  [[nodiscard]] std::optional<size_t>
  num_matrices(const std::filesystem::path &input_path) const;

  // files: all input files, in order, with their number of matrices
  void set_new_block_mapping(
    const std::vector<std::pair<std::filesystem::path, size_t>> &files);

  [[nodiscard]] size_t num_reused_blocks() const;

  // Old pmp_info entry for old block old_index, renumbered to new_index
  [[nodiscard]] std::string
  renumbered_pmp_info_entry(size_t old_index, size_t new_index) const;

  [[nodiscard]] std::filesystem::path normalization_path() const;
};
