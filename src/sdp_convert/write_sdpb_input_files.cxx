#include "Dual_Constraint_Group.hxx"
#include "../set_stream_precision.hxx"

#include <boost/filesystem.hpp>

void write_control(const boost::filesystem::path &output_dir,
                   const size_t &num_blocks,
                   const std::vector<std::string> &command_arguments);

void write_objectives(const boost::filesystem::path &output_dir,
                      const El::BigFloat &objective_const,
                      const std::vector<El::BigFloat> &dual_objective_b);

void write_bilinear_bases(boost::filesystem::ofstream &output_stream,
                          const Dual_Constraint_Group &group);

void write_blocks(boost::filesystem::ofstream &output_stream,
                  const Dual_Constraint_Group &group);

void write_primal_objective_c(boost::filesystem::ofstream &output_stream,
                              const Dual_Constraint_Group &group);

void write_free_var_matrix(boost::filesystem::ofstream &output_stream,
                           const size_t &dual_objectives_b_size,
                           const Dual_Constraint_Group &group);

void write_sdpb_input_files(
  const boost::filesystem::path &output_dir, const int &rank,
  const size_t &num_blocks, const std::vector<std::string> &command_arguments,
  const El::BigFloat &objective_const,
  const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  boost::filesystem::create_directories(output_dir);
  if(rank == 0)
    {
      write_control(output_dir, num_blocks, command_arguments);
      write_objectives(output_dir, objective_const, dual_objective_b);
    }
  for(auto &group : dual_constraint_groups)
    {
      const boost::filesystem::path output_path(
        output_dir
        / ("block_" + std::to_string(group.block_index) + ".json"));
      boost::filesystem::ofstream output_stream(output_path);
      set_stream_precision(output_stream);
      output_stream << "{\n";

      write_blocks(output_stream, group);
      write_bilinear_bases(output_stream, group);
      write_primal_objective_c(output_stream, group);
      write_free_var_matrix(output_stream, dual_objective_b.size(), group);
      output_stream << "}\n";
      if(!output_stream.good())
        {
          throw std::runtime_error("Error when writing to: "
                                   + output_path.string());
        }
    }
}
