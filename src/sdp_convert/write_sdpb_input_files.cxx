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

void print_matrix_sizes(
  const int &rank, const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups);

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
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups,
  const bool debug)
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
  if(debug)
    {
      print_matrix_sizes(rank, dual_objective_b, dual_constraint_groups);
    }
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

void print_matrix_sizes(
  const int &rank, const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups)
{
  if(rank == 0)
    {
      El::Output("---------------------");
      El::Output("Matrix sizes and RAM estimates:");
      El::Output("---------------------");
    }

  size_t N = dual_objective_b.size();
  size_t my_P = 0; // P' for a given group

  // Hereafter XXX_elements means number of nonzero elements in XXX
  size_t my_constraint_matrix_elements = 0;
  size_t my_bilinear_bases_elements = 0;
  size_t my_bilinear_pairing_block_elements = 0;
  size_t my_psd_elements = 0;
  size_t my_schur_elements = 0;
  for(auto &group : dual_constraint_groups)
    {
      my_P += group.constraint_constants.size();

      my_constraint_matrix_elements += group.constraint_matrix.MemorySize();
      my_bilinear_bases_elements += group.bilinear_bases[0].MemorySize()
                                    + group.bilinear_bases[1].MemorySize();

      // variables stored in Block_Info
      size_t dimensions = group.dim;
      size_t num_points = group.degree + 1; // see write_blocks()

      // See Block_Info::schur_block_sizes()
      size_t schur_width = num_points * dimensions * (dimensions + 1) / 2;
      my_schur_elements += schur_width * schur_width;

      // See Block_Info::bilinear_pairing_block_sizes()
      // two blocks, each one has same size
      size_t bilinear_pairing_block_width = num_points * dimensions;
      my_bilinear_pairing_block_elements
        += 2 * bilinear_pairing_block_width * bilinear_pairing_block_width;

      // See Block_Info::psd_matrix_block_sizes()
      // for each group we have two psd square matrix blocks, with width
      // psd_even and psd_odd
      auto psd_even = dimensions * ((num_points + 1) / 2);
      auto psd_odd = dimensions * num_points - psd_even;
      my_psd_elements += psd_even * psd_even + psd_odd * psd_odd;
    }
  if(my_P * N != my_constraint_matrix_elements)
    {
      El::LogicError("sum(P'*N) != #(B bands): P'=", my_P, ", N=", N,
                     "#(B bands)=", my_constraint_matrix_elements);
    }

  // NB: Reduce should be called on all ranks, its result used only on rank 0!
  auto reduce_sum = [](size_t item) {
    return El::mpi::Reduce(item, El::mpi::SUM, 0, El::mpi::COMM_WORLD);
  };

  size_t B_matrix_elements = reduce_sum(my_constraint_matrix_elements);
  size_t bilinear_pairing_block_elements
    = reduce_sum(my_bilinear_pairing_block_elements);
  size_t psd_blocks_elements = reduce_sum(my_psd_elements);
  size_t schur_elements = reduce_sum(my_schur_elements);
  size_t P = reduce_sum(my_P);
  size_t bilinear_bases_elements = reduce_sum(my_bilinear_bases_elements);

  if(rank == 0)
    {
      size_t big_float_bytes = El::BigFloat(1.1).SerializedSize();

      size_t Q_matrix_size = N * N;
      size_t total_no_Q = 2 * B_matrix_elements + 5 * psd_blocks_elements
                          + 2 * schur_elements
                          + 2 * bilinear_pairing_block_elements;

      std::vector<std::pair<std::string, size_t>> sizes{
        {"P (primal objective)", P},
        {"N (dual objective)", N},
        {"B matrix (PxN)", B_matrix_elements},
        {"Q matrix (NxN)", Q_matrix_size},
        {"Bilinear bases", bilinear_bases_elements},
        {"Bilinear pairing blocks", bilinear_pairing_block_elements},
        {"PSD blocks", psd_blocks_elements},
        {"Schur (PxP block diagonal)", schur_elements},
        {"Total (no Q) = 2#(B) + 5#(PSD) + 2#(S) + 2#(Bilinear pairing)",
         total_no_Q},
      };

      El::Output("BigFloat, bytes: ", big_float_bytes);
      for(auto &[key, value] : sizes)
        {
          El::Output(key, ", elements: ", value);
          El::Output(key, ", bytes: ", value * big_float_bytes);
        }
      double GB_per_element = (double)big_float_bytes / 1024 / 1024 / 1024;
      El::Output("Total RAM (no Q), GB: ", total_no_Q * GB_per_element);
      El::Output("Q, GB: ", Q_matrix_size * GB_per_element);
      El::Output("NB: Q is copied over each core group, i.e. "
                 "(nodes*procsPerNode/procGranularity) times");
      El::Output("---------------------");
    }
}
