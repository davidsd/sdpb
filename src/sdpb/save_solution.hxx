#pragma once

#include "SDPB_Parameters.hxx"
#include "pmp/max_normalization_index.hxx"
#include "sdp_solve/sdp_solve.hxx"
#include "sdpa_solve/sdpa_solve.hxx"
#include "sdpb_util/ostream/set_stream_precision.hxx"
#include "sdpb_util/write_matrix.hxx"

inline void write_block_vector(const Block_Vector &v,
                               const std::vector<size_t> &block_indices,
                               const std::filesystem::path &prefix)
{
  for(size_t b = 0; b < v.blocks.size(); ++b)
    {
      auto &block = v.blocks[b];
      const auto block_index = block_indices.at(b);
      const auto path = prefix.string() + El::BuildString(block_index, ".txt");
      write_distmatrix(block, path);
    }
}

inline void write_block_matrix(const Block_Diagonal_Matrix &m,
                               const std::vector<size_t> &block_indices,
                               const std::filesystem::path &prefix)
{
  for(size_t b = 0; b < m.blocks.size(); ++b)
    {
      auto &block = m.blocks[b];
      const auto block_index = block_indices.at(b);
      const auto path = prefix.string() + El::BuildString(block_index, ".txt");
      write_distmatrix(block, path);
    }
}
inline void write_paired_block_matrix(const Paired_Block_Diagonal_Matrix &m,
                                      const std::vector<size_t> &block_indices,
                                      const std::filesystem::path &prefix)
{
  for(size_t b = 0; b < m.blocks.size(); ++b)
    {
      auto &block = m.blocks[b];
      if(block.Height() == 0)
        continue;
      const auto block_index = block_indices.at(b / 2);
      const auto parity = b % 2;
      const auto path
        = prefix.string() + El::BuildString(2 * block_index + parity, ".txt");
      write_distmatrix(block, path);
    }
}

inline El::BigFloat primal_error(const SDP_Solver &solver)
{
  return solver.primal_error();
}
inline El::BigFloat primal_error(const Sdpb::Sdpa::SDP_Solver &solver)
{
  return solver.primal_error;
}

template <class TSolver>
void write_out_txt(const TSolver &solver,
                   const SDP_Solver_Terminate_Reason &terminate_reason,
                   const int64_t &solver_runtime,
                   const std::filesystem::path &out_directory,
                   const Verbosity &verbosity)
{
  // Internally, El::Print() sync's everything to the root core and
  // outputs it from there.  So do not actually open the file on
  // anything but the root node.

  std::ofstream out_stream;
  if(El::mpi::Rank() == 0)
    {
      if(verbosity >= Verbosity::regular)
        {
          std::cout << "Saving solution to      : " << out_directory << '\n';
        }
      create_directories(out_directory);
      const auto output_path = out_directory / "out.txt";
      out_stream.open(output_path);
      set_stream_precision(out_stream);
      out_stream << "terminateReason = \"" << terminate_reason << "\";\n"
                 << "primalObjective = " << solver.primal_objective << ";\n"
                 << "dualObjective   = " << solver.dual_objective << ";\n"
                 << "dualityGap      = " << solver.duality_gap << ";\n"
                 << "primalError     = " << primal_error(solver) << ";\n"
                 << "dualError       = " << solver.dual_error << ";\n"
                 << "Solver runtime  = " << solver_runtime << ";\n";
      ASSERT(out_stream.good(), "Error when writing to: ", output_path);
    }
}

[[nodiscard]]
inline El::Matrix<El::BigFloat>
get_z(const El::Matrix<El::BigFloat> &y,
      const std::vector<El::BigFloat> &normalization)
{
  El::Matrix<El::BigFloat> z(y.Height() + 1, y.Width());
  ASSERT_EQUAL(normalization.size(), y.Height() + 1);

  auto max_index = max_normalization_index(normalization);
  ASSERT(normalization.at(max_index) != El::BigFloat(0));
  // To construct z, we take y
  // and insert a new element into it at max_index position.
  // The value of this element is chosen to satisfy
  // the normalization condition n.z == 1.
  // NB: this procedure (in particular, choice of max_index)
  // is the inverse of conversion
  // Polynomial_Matrix_Program.objective -> Output_SDP.dual_objective_b
  // in Outout_SDP constructor
  for(int i = 0; i < max_index; ++i)
    z(i, 0) = y(i, 0);
  z(max_index, 0) = 0; // to be updated below
  for(int i = max_index; i < y.Height(); ++i)
    z(i + 1, 0) = y(i, 0);

  // Calculate n.z
  // we convert sdt::vector to El::Matrix to use it in El::Dot
  El::Matrix<El::BigFloat> normalization_view;
  normalization_view.LockedAttach(normalization.size(), 1,
                                  normalization.data(), normalization.size());
  // n.z - n[max_index]*z[max_index]
  auto nz = El::Dot(normalization_view, z);
  z(max_index, 0) = (1 - nz) / normalization.at(max_index);
  // now n.z == 1
  return z;
}

inline void
save_solution(const SDP_Solver &solver, const Block_Info &block_info,
              const SDP &sdp, const SDPB_Parameters &parameters,
              const SDP_Solver_Terminate_Reason &terminate_reason,
              const int64_t &solver_runtime)
{
  const auto &out_directory = parameters.out_directory;
  const auto &write_solution = parameters.write_solution;
  const auto &block_indices = block_info.block_indices;

  write_out_txt(solver, terminate_reason, solver_runtime, out_directory,
                parameters.verbosity);
  if(write_solution.vector_y || write_solution.vector_z)
    {
      const auto y_path = out_directory / "y.txt";
      const auto z_path = out_directory / "z.txt";
      // y is duplicated among blocks, so only need to print out copy
      // from the first block of rank 0.
      const auto &y_dist = solver.y.blocks.at(0);
      if(y_dist.Root() == 0)
        {
          // Copy from all ranks owning a block to rank zero
          El::DistMatrix<El::BigFloat, El::CIRC, El::CIRC> y_circ(y_dist);
          ASSERT_EQUAL(y_circ.Root(), 0);
          if(El::mpi::Rank() == 0)
            {
              // local matrix
              const El::Matrix<El::BigFloat> &y = y_circ.Matrix();
              ASSERT_EQUAL(y.Height(), y_dist.Height());
              ASSERT_EQUAL(y.Width(), 1);
              if(write_solution.vector_y)
                {
                  write_matrix(y, y_path);
                }
              if(write_solution.vector_z)
                {
                  ASSERT(sdp.normalization.has_value());
                  const auto z = get_z(y, sdp.normalization.value());
                  write_matrix(z, z_path);
                }
            }
        }
      if(El::mpi::Rank() == 0)
        {
          // Sanity check.
          // Assertion will fail if something went completely wrong
          // and (y_dist.Root() == 0) is false for all blocks
          if(write_solution.vector_y)
            ASSERT(exists(y_path), y_path);
          if(write_solution.vector_z)
            ASSERT(exists(z_path), z_path);
        }
    }

  if(write_solution.vector_x)
    write_block_vector(solver.x, block_indices, out_directory / "x_");
  if(write_solution.matrix_X)
    write_paired_block_matrix(solver.X, block_indices,
                              out_directory / "X_matrix_");
  if(write_solution.matrix_Y)
    write_paired_block_matrix(solver.Y, block_indices,
                              out_directory / "Y_matrix_");
}

inline void save_solution(const Sdpb::Sdpa::SDP_Solver &solver,
                          const Sdpb::Sdpa::Block_Info &block_info,
                          [[maybe_unused]] const Sdpb::Sdpa::SDP &sdp,
                          const SDPB_Parameters &parameters,
                          const SDP_Solver_Terminate_Reason &terminate_reason,
                          const int64_t &solver_runtime)
{
  const auto &out_directory = parameters.out_directory;
  const auto &write_solution = parameters.write_solution;
  const auto &block_indices = block_info.block_indices;

  write_out_txt(solver, terminate_reason, solver_runtime, out_directory,
                parameters.verbosity);
  if(El::mpi::Rank() == 0)
    {
      // TODO Use another Write_Solution type?
      if(write_solution.vector_y)
        PRINT_WARNING("--writeSolution=y is not available for SDPA problems");
      if(write_solution.vector_z)
        PRINT_WARNING("--writeSolution=z is not available for SDPA problems");
    }

  if(write_solution.vector_x)
    write_distmatrix(solver.x, out_directory / "x.txt");
  if(write_solution.matrix_X)
    write_block_matrix(solver.X, block_indices, out_directory / "X_matrix_");
  if(write_solution.matrix_Y)
    write_block_matrix(solver.Y, block_indices, out_directory / "Y_matrix_");
}
