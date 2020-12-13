#include "Dual_Constraint_Group.hxx"

#include <boost/filesystem.hpp>

#include <vector>

void write_control(const boost::filesystem::path &output_dir,
                   const int &num_procs,
                   const std::vector<std::string> &command_arguments);

void write_objectives(const boost::filesystem::path &output_dir,
                      const El::BigFloat &objective_const,
                      const std::vector<El::BigFloat> &dual_objective_b);

void write_bilinear_bases(
  const boost::filesystem::path &output_dir, const int &rank,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups);

void write_blocks(
  const boost::filesystem::path &output_dir, const int &rank,
  const std::vector<size_t> &indices,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups);

void write_primal_objective_c(
  const boost::filesystem::path &output_dir,
  const std::vector<size_t> &indices,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups);

void write_free_var_matrix(
  const boost::filesystem::path &output_dir,
  const std::vector<size_t> &indices, const size_t &dual_objectives_b_size,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups);

void write_sdpb_input_files(
  const boost::filesystem::path &output_dir, const int &rank,
  const int &num_procs, const std::vector<std::string> &command_arguments,
  const std::vector<size_t> &indices, const El::BigFloat &objective_const,
  const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  boost::filesystem::create_directories(output_dir);
  if(rank == 0)
    {
      write_control(output_dir, num_procs, command_arguments);
      write_objectives(output_dir, objective_const, dual_objective_b);
    }
  write_bilinear_bases(output_dir, rank, dual_constraint_groups);
  write_blocks(output_dir, rank, indices, dual_constraint_groups);
  write_primal_objective_c(output_dir, indices, dual_constraint_groups);
  write_free_var_matrix(output_dir, indices, dual_objective_b.size(),
                        dual_constraint_groups);
}
