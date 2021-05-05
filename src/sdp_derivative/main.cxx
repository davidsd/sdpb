#include "DSDP_Parameters.hxx"
#include "../sdp_solve.hxx"

#include "../set_stream_precision.hxx"

#include <El.hpp>

void compute_B_pseudoinverse(const std::vector<size_t> &block_offsets,
                             const std::vector<size_t> &block_indices,
                             const El::Grid &grid, const Block_Matrix &B,
                             Block_Matrix &B_pseudoinverse);
void Axpy(const El::BigFloat &alpha, const SDP &new_sdp, SDP &delta_sdp);

void load_solution(const boost::filesystem::path &checkpoint_path,
                   Block_Vector &x, Block_Vector &y);
El::BigFloat
compute_approximate_objective(const SDP &sdp, const SDP &d_sdp,
                              const Block_Vector &x, const Block_Vector &y,
                              const Block_Matrix &B_pseudoinverse);

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
                            parameters.proc_granularity, parameters.verbosity);

      El::Grid grid(block_info.mpi_comm.value);
      std::cout << "param: "
                << parameters << "\n";
      
      set_stream_precision(std::cout);
      SDP sdp(parameters.sdp_path, block_info, grid);

      std::vector<size_t> block_offsets(sdp.free_var_matrix.blocks.size() + 1,
                                        0);
      for(size_t p(0); p < sdp.free_var_matrix.blocks.size(); ++p)
        {
          block_offsets[p + 1]
            = block_offsets[p] + sdp.free_var_matrix.blocks[p].Height();
        }
      Block_Matrix B_pseudoinverse(sdp.free_var_matrix);
      compute_B_pseudoinverse(block_offsets, block_info.block_indices, grid,
                              sdp.free_var_matrix, B_pseudoinverse);

      SDP new_sdp(parameters.new_sdp_path, block_info, grid), d_sdp(new_sdp);
      Axpy(El::BigFloat(-1), sdp, d_sdp);

      // El::Print(sdp.free_var_matrix.blocks.front(),"B");
      // // El::Print(new_sdp.free_var_matrix.blocks.front(),"\nB_new");
      // El::Print(B_pseudoinverse.blocks.front(),"\nB_pinv");
      // El::Print(d_sdp.free_var_matrix.blocks.front(),"\ndB");
      // // El::Print(d_sdp.primal_objective_c.blocks.front(),"\ndc");
      // // El::Print(d_sdp.dual_objective_b,"\ndb");
      // std::cout << "\n";

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

      El::BigFloat new_objective(
        compute_approximate_objective(sdp, d_sdp, x, y, B_pseudoinverse));
      std::cout << "objective: "
                << new_objective
                << "\n";
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
