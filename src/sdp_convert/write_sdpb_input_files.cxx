#include "Dual_Constraint_Group.hxx"
#include "../set_stream_precision.hxx"

#include <archive_writer.hpp>
#include <archive_exception.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

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
  const boost::filesystem::path &output_path, const int &rank,
  const size_t &num_blocks, const std::vector<std::string> &command_arguments,
  const El::BigFloat &objective_const,
  const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  boost::filesystem::path temp_dir(output_path);
  temp_dir+="_temp";
  boost::filesystem::create_directories(temp_dir);
  if(rank == 0)
    {
      write_control(temp_dir, num_blocks, command_arguments);
      write_objectives(temp_dir, objective_const, dual_objective_b);
    }
  for(auto &group : dual_constraint_groups)
    {
      const boost::filesystem::path block_path(
        temp_dir / ("block_" + std::to_string(group.block_index) + ".json"));
      boost::filesystem::ofstream output_stream(block_path);
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
                                   + block_path.string());
        }
    }
  El::mpi::Barrier(El::mpi::COMM_WORLD);
  if(rank == 0)
    {
      boost::filesystem::ofstream zip_stream(output_path);
      ns_archive::writer writer(
        ns_archive::writer::make_writer<ns_archive::ns_writer::format::_ZIP>(
          zip_stream, 10240));

      const boost::filesystem::directory_iterator end;
      for(auto &entry : boost::filesystem::directory_iterator(temp_dir))
        {
          {
            boost::filesystem::ifstream stream(entry.path());
            ns_archive::entry out_entry(stream);
            out_entry.set_header_value_pathname(
                                                entry.path().filename().string());
            writer.add_entry(out_entry);
          }
          boost::filesystem::remove(entry.path());
        }
    }
  boost::filesystem::remove(temp_dir);
}
