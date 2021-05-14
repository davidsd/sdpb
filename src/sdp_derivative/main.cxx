#include "DSDP_Parameters.hxx"
#include "../sdp_solve.hxx"

#include "../set_stream_precision.hxx"

#include <El.hpp>

void Axpy(const El::BigFloat &alpha, const SDP &new_sdp, SDP &delta_sdp);

void load_solution(const boost::filesystem::path &checkpoint_path,
                   Block_Vector &x, Block_Vector &y);
El::BigFloat
compute_approximate_objective(const Block_Info &block_info,
                              const El::Grid &grid, const SDP &sdp,
                              const SDP &d_sdp, const Block_Vector &x,
                              const Block_Vector &y,
                              const Block_Diagonal_Matrix &X,
                              const Block_Diagonal_Matrix &Y);

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  try
    {
      DSDP_Parameters parameters(argc, argv);
      if(!parameters.is_valid())
        {
          return 0;
        }
      El::gmp::SetPrecision(parameters.precision);
      Block_Info block_info(parameters.sdp_path, parameters.solution_dir,
                            parameters.procs_per_node,
                            parameters.proc_granularity, Verbosity::none);

      El::Grid grid(block_info.mpi_comm.value);
      set_stream_precision(std::cout);
      SDP sdp(parameters.sdp_path, block_info, grid);

      std::vector<size_t> block_offsets(sdp.free_var_matrix.blocks.size() + 1,
                                        0);
      for(size_t p(0); p < sdp.free_var_matrix.blocks.size(); ++p)
        {
          block_offsets[p + 1]
            = block_offsets[p] + sdp.free_var_matrix.blocks[p].Height();
        }

      SDP new_sdp(parameters.new_sdp_path, block_info, grid), d_sdp(new_sdp);
      Axpy(El::BigFloat(-1), sdp, d_sdp);

      Block_Vector x(block_info.schur_block_sizes(), block_info.block_indices,
                     block_info.num_points.size(), grid),
        y(std::vector<size_t>(block_info.num_points.size(),
                              sdp.dual_objective_b.Height()),
          block_info.block_indices, block_info.num_points.size(), grid);
      Block_Diagonal_Matrix X(block_info.psd_matrix_block_sizes(),
                              block_info.block_indices,
                              block_info.num_points.size(), grid), Y(X);
      for(size_t block = 0; block != block_info.block_indices.size(); ++block)
        {
          size_t block_index(block_info.block_indices.at(block));
          read_text_block(x.blocks.at(block), parameters.solution_dir, "x_",
                          block_index);
          read_text_block(y.blocks.at(block),
                          parameters.solution_dir / "y.txt");
          for(size_t psd_block(0); psd_block < 2; ++psd_block)
            {
              // Constant constraints have empty odd parity blocks, so we do
              // not need to load them.
              if(X.blocks.at(2 * block + psd_block).Height() != 0)
                {
                  const size_t psd_index(2 * block_index + psd_block);
                  read_text_block(X.blocks.at(2 * block + psd_block),
                                  parameters.solution_dir, "X_matrix_",
                                  psd_index);
                  read_text_block(Y.blocks.at(2 * block + psd_block),
                                  parameters.solution_dir, "Y_matrix_",
                                  psd_index);
                }
            }
        }

      El::BigFloat new_objective(
                                 compute_approximate_objective(block_info, grid, sdp, d_sdp, x, y, X, Y));
      std::cout << "{\"approx_objective\": \"" << new_objective << "\"}\n";
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
