#include "Dual_Constraint_Group.hxx"
#include "../set_stream_precision.hxx"

#include "Archive_Writer.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

void write_control(const boost::filesystem::path &output_dir,
                   const size_t &num_blocks,
                   const std::vector<std::string> &command_arguments);

void write_objectives(const boost::filesystem::path &output_dir,
                      const El::BigFloat &objective_const,
                      const std::vector<El::BigFloat> &dual_objective_b);

void write_bilinear_bases(std::ostream &output_stream,
                          const Dual_Constraint_Group &group);

void write_blocks(std::ostream &output_stream,
                  const Dual_Constraint_Group &group);

void write_primal_objective_c(std::ostream &output_stream,
                              const Dual_Constraint_Group &group);

void write_free_var_matrix(std::ostream &output_stream,
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
  temp_dir += "_temp";
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

      // Use gzip with no compression to get a CRC
      boost::iostreams::filtering_ostream output_stream;
      output_stream.push(
        boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(0)));
      output_stream.push(boost::iostreams::file_sink(block_path.string()));
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
      std::vector<boost::filesystem::path> paths_to_remove;
      Archive_Writer writer(output_path);
      for(auto &directory_entry : boost::filesystem::directory_iterator(temp_dir))
        {
          paths_to_remove.emplace_back(directory_entry.path());
          auto &path(paths_to_remove.back());
          boost::iostreams::filtering_stream<boost::iostreams::input> input_stream;
          input_stream.push(boost::iostreams::gzip_decompressor());
          input_stream.push(boost::iostreams::file_source(path.string()));
  
          writer.write_entry(Archive_Entry(path), input_stream);
        }
      paths_to_remove.emplace_back(temp_dir);
      for(auto &path : paths_to_remove)
        {
          boost::filesystem::remove(path);
        }
    }
}
