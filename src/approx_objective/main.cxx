#include "Approx_Parameters.hxx"
#include "../sdp_solve.hxx"

#include "../set_stream_precision.hxx"

#include <El.hpp>

void setup_solver(const Block_Info &block_info, const El::Grid &grid,
                  const SDP &sdp, const boost::filesystem::path &solution_dir,
                  Block_Diagonal_Matrix &schur_complement_cholesky,
                  Block_Matrix &schur_off_diagonal,
                  El::DistMatrix<El::BigFloat> &Q);

void write_solver_state(const std::vector<size_t> &block_indices,
                        const boost::filesystem::path &solution_dir,
                        const Block_Diagonal_Matrix &schur_complement_cholesky,
                        const Block_Matrix &schur_off_diagonal,
                        const El::DistMatrix<El::BigFloat> &Q);

std::vector<std::pair<std::string, El::BigFloat>>
compute_approximate_objectives(
  const Block_Info &block_info, const El::Grid &grid, const SDP &sdp,
  const Block_Vector &x, const Block_Vector &y,
  const Block_Diagonal_Matrix &schur_complement_cholesky,
  const Block_Matrix &schur_off_diagonal,
  const El::DistMatrix<El::BigFloat> &Q,
  const boost::filesystem::path &input_path);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      Approx_Parameters parameters(argc, argv);
      if(!parameters.is_valid())
        {
          return 0;
        }
      El::gmp::SetPrecision(parameters.precision);
      Block_Info block_info(parameters.sdp_path, parameters.solution_dir,
                            parameters.procs_per_node,
                            parameters.proc_granularity, Verbosity::none);

      El::Grid grid(block_info.mpi_comm.value);
      SDP sdp(parameters.sdp_path, block_info, grid);

      std::vector<size_t> block_offsets(sdp.free_var_matrix.blocks.size() + 1,
                                        0);
      for(size_t p(0); p < sdp.free_var_matrix.blocks.size(); ++p)
        {
          block_offsets[p + 1]
            = block_offsets[p] + sdp.free_var_matrix.blocks[p].Height();
        }

      Block_Vector x(block_info.schur_block_sizes(), block_info.block_indices,
                     block_info.num_points.size(), grid),
        y(std::vector<size_t>(block_info.num_points.size(),
                              sdp.dual_objective_b.Height()),
          block_info.block_indices, block_info.num_points.size(), grid);
      for(size_t block = 0; block != block_info.block_indices.size(); ++block)
        {
          size_t block_index(block_info.block_indices.at(block));
          read_text_block(x.blocks.at(block), parameters.solution_dir, "x_",
                          block_index);
          read_text_block(y.blocks.at(block),
                          parameters.solution_dir / "y.txt");
        }

      Block_Diagonal_Matrix schur_complement_cholesky(
        block_info.schur_block_sizes(), block_info.block_indices,
        block_info.num_points.size(), grid);
      Block_Matrix schur_off_diagonal(sdp.free_var_matrix);
      El::DistMatrix<El::BigFloat> Q(sdp.dual_objective_b.Height(),
                                     sdp.dual_objective_b.Height());

      setup_solver(block_info, grid, sdp, parameters.solution_dir,
                   schur_complement_cholesky, schur_off_diagonal, Q);
      if(parameters.write_solver_state)
        {
          write_solver_state(block_info.block_indices, parameters.solution_dir,
                             schur_complement_cholesky, schur_off_diagonal, Q);
        }
      std::vector<std::pair<std::string, El::BigFloat>> approx_objectives;
      if(!parameters.new_sdp_path.empty())
        {
          approx_objectives = compute_approximate_objectives(
            block_info, grid, sdp, x, y, schur_complement_cholesky,
            schur_off_diagonal, Q, parameters.new_sdp_path);
        }
      if(El::mpi::Rank() == 0)
        {
          set_stream_precision(std::cout);
          std::cout << "[\n";
          for(auto iter(approx_objectives.begin());
              iter != approx_objectives.end(); ++iter)
            {
              if(iter != approx_objectives.begin())
                {
                  std::cout << ",\n";
                }
              std::cout << "  {\n"
                        << "    \"path\": \"" << iter->first << "\",\n"
                        << "    \"objective\": \"" << iter->second << "\"\n"
                        << "  }";
            }
          std::cout << "\n]\n";
        }
    }
  catch(std::exception &e)
    {
      El::ReportException(e);
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
  catch(...)
    {
      El::mpi::Abort(El::mpi::COMM_WORLD, 1);
    }
}
