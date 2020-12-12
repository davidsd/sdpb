#include "Dual_Constraint_Group.hxx"
#include "write_vector.hxx"

void write_blocks(
  const boost::filesystem::path &output_dir, const int &rank,
  const int &num_procs, const std::vector<size_t> &indices,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  std::vector<size_t> dimensions, degrees, schur_block_sizes,
    psd_matrix_block_sizes, bilinear_pairing_block_sizes;
  for(auto &g : dual_constraint_groups)
    {
      dimensions.push_back(g.dim);
      degrees.push_back(g.degree);

      schur_block_sizes.push_back((g.dim * (g.dim + 1) / 2) * (g.degree + 1));

      // sdp.bilinear_bases is the concatenation of the g.bilinear_bases.
      // The matrix Y is a BlockDiagonalMatrix built from the
      // concatenation of the blocks for each individual
      // Dual_Constraint_Group.  sdp.blocks[j] = {b1, b2, ... } contains
      // the indices for the blocks of Y corresponding to the j-th
      // group.
      for(auto &b : g.bilinear_bases)
        {
          psd_matrix_block_sizes.push_back(b.Height() * g.dim);
          bilinear_pairing_block_sizes.push_back(b.Width() * g.dim);
        }
    }

  const boost::filesystem::path output_path(
    output_dir / ("blocks_" + std::to_string(rank) + ".json"));
  boost::filesystem::ofstream output_stream(output_path);
  output_stream << "{\n  \"num_procs\": " << num_procs << ",\n  ";
  write_vector(output_stream, indices, "indices");
  output_stream << ",\n  ";
  write_vector(output_stream, dimensions, "dimensions");
  output_stream << ",\n  ";
  write_vector(output_stream, degrees, "degrees");
  output_stream << ",\n  ";
  write_vector(output_stream, schur_block_sizes, "schur_block_sizes");
  output_stream << ",\n  ";
  write_vector(output_stream, psd_matrix_block_sizes,
               "psd_matrix_block_sizes");
  output_stream << ",\n  ";
  write_vector(output_stream, bilinear_pairing_block_sizes,
               "bilinear_pairing_block_sizes");
  output_stream << "\n}\n";
  if(!output_stream.good())
    {
      throw std::runtime_error("Error when writing to: "
                               + output_path.string());
    }
}
