#include "Dual_Constraint_Group.hxx"
#include "byte_counter.hxx"
#include "Archive_Writer.hxx"
#include "../set_stream_precision.hxx"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

size_t write_control(const boost::filesystem::path &output_dir,
                     const size_t &num_blocks,
                     const std::vector<std::string> &command_arguments);

size_t write_objectives(const boost::filesystem::path &output_dir,
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

namespace
{
  void archive_gzipped_file(const boost::filesystem::path &path,
                            const int64_t &num_bytes, Archive_Writer &writer)
  {
    boost::iostreams::filtering_stream<boost::iostreams::input> input_stream;
    input_stream.push(boost::iostreams::gzip_decompressor());
    input_stream.push(boost::iostreams::file_source(path.string()));

    writer.write_entry(Archive_Entry(path, num_bytes), input_stream);
  }
}

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
  size_t num_control_bytes(0), num_objectives_bytes(0);
  if(rank == 0)
    {
      num_control_bytes
        = write_control(temp_dir, num_blocks, command_arguments);
      num_objectives_bytes
        = write_objectives(temp_dir, objective_const, dual_objective_b);
    }
  // We use size_t rather than std::streamsize because MPI treats
  // std::streamsize as an MPI_LONG_INT and then can not MPI_Reduce
  // over it.
  std::vector<size_t> block_file_sizes(num_blocks, 0);
  for(auto &group : dual_constraint_groups)
    {
      const boost::filesystem::path block_path(
        temp_dir / ("block_" + std::to_string(group.block_index) + ".json"));

      byte_counter counter;
      {
        boost::iostreams::filtering_ostream output_stream;
        output_stream.push(boost::ref(counter));
        // Use gzip with no compression to get a CRC
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
      block_file_sizes.at(group.block_index) = counter.num_bytes;
    }
  El::mpi::Reduce(block_file_sizes.data(), block_file_sizes.size(),
                  El::mpi::SUM, 0, El::mpi::COMM_WORLD);
  if(rank == 0)
    {
      Archive_Writer writer(output_path);
      const boost::filesystem::path control_path(temp_dir / "control.json"),
        objectives_path(temp_dir / "objectives.json");
      archive_gzipped_file(control_path, num_control_bytes, writer);
      boost::filesystem::remove(control_path);

      archive_gzipped_file(objectives_path, num_objectives_bytes, writer);
      boost::filesystem::remove(objectives_path);

      for(size_t block_index(0); block_index != block_file_sizes.size();
          ++block_index)
        {
          const boost::filesystem::path block_path(
            temp_dir / ("block_" + std::to_string(block_index) + ".json"));
          archive_gzipped_file(block_path, block_file_sizes.at(block_index),
                               writer);
          boost::filesystem::remove(block_path);
        }
      boost::filesystem::remove(temp_dir);
    }
}
