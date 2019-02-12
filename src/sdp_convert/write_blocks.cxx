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
      if(g.bilinear_bases.size() != 2)
        {
          throw std::runtime_error("Wrong number of elements in "
                                   "dualConstraintGroups::bilinear_bases.  "
                                   "Expected 2 but found: "
                                   + std::to_string(g.bilinear_bases.size()));
        }
      for(auto &b : g.bilinear_bases)
        {
          psd_matrix_block_sizes.push_back(b.Height() * g.dim);
          bilinear_pairing_block_sizes.push_back(b.Width() * g.dim);
        }
    }

  boost::filesystem::ofstream output_stream(
    output_dir / ("blocks." + std::to_string(rank)));
  output_stream << num_procs << "\n";
  write_vector(output_stream, indices);
  write_vector(output_stream, dimensions);
  write_vector(output_stream, degrees);
  write_vector(output_stream, schur_block_sizes);
  write_vector(output_stream, psd_matrix_block_sizes);
  write_vector(output_stream, bilinear_pairing_block_sizes);
}
