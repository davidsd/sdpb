#include "Dual_Constraint_Group.hxx"
#include "byte_counter.hxx"
#include "Archive_Writer.hxx"
#include "write_sdp.hxx"
#include "Output_SDP/Output_SDP.hxx"
#include "sdpb_util/assert.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <filesystem>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace fs = std::filesystem;

size_t
write_control_json(std::ostream &output_stream, const size_t &num_blocks,
                   const std::vector<std::string> &command_arguments);

void write_objectives_json(std::ostream &output_stream,
                           const El::BigFloat &objective_const,
                           const std::vector<El::BigFloat> &dual_objective_b);

void write_block_info_json(std::ostream &output_stream,
                           const Dual_Constraint_Group &group);

void write_block_data(std::ostream &output_stream,
                      const Dual_Constraint_Group &group,
                      Block_File_Format format);

void print_matrix_sizes(
  const int &rank, const std::vector<El::BigFloat> &dual_objective_b,
  const std::vector<Dual_Constraint_Group> &dual_constraint_groups);

namespace
{
  size_t write_data_and_count_bytes(
    const fs::path &output_path,
    const std::function<void(std::ostream &)> &write_data, const bool zip,
    const bool binary = false)
  {
    byte_counter counter;
    boost::iostreams::filtering_ostream output_stream;
    output_stream.push(boost::ref(counter));
    if(zip)
      {
        // Use gzip with no compression to get a CRC
        output_stream.push(
          boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(0)));
      }
    auto openmode = std::ios::out;
    if(binary)
      openmode |= std::ios::binary;
    output_stream.push(boost::iostreams::file_sink(output_path.string()),
                       openmode);
    if(!binary)
      set_stream_precision(output_stream);
    write_data(output_stream);
    ASSERT(output_stream.good(), "Error when writing to: ", output_path);
    return counter.num_bytes;
  }

  fs::path get_block_info_path(const fs::path &temp_dir, size_t block_index)
  {
    return temp_dir / El::BuildString("block_info_", block_index, ".json");
  }

  fs::path get_block_data_path(const fs::path &temp_dir, size_t block_index,
                               Block_File_Format format)
  {
    std::string name = El::BuildString("block_data_", block_index);
    switch(format)
      {
      case json: name += ".json"; break;
      case bin: name += ".bin"; break;
      default: RUNTIME_ERROR("Unsupported Block_File_Format: ", format);
      }
    return temp_dir / name;
  }

  fs::path get_control_path(const fs::path &temp_dir)
  {
    return temp_dir / "control.json";
  }
  fs::path get_objectives_path(const fs::path &temp_dir)
  {
    return temp_dir / "objectives.json";
  }

  void archive_gzipped_file(const fs::path &path, const int64_t &num_bytes,
                            Archive_Writer &writer)
  {
    boost::iostreams::filtering_stream<boost::iostreams::input> input_stream;
    input_stream.push(boost::iostreams::gzip_decompressor());
    input_stream.push(boost::iostreams::file_source(path.string()));

    writer.write_entry(Archive_Entry(path, num_bytes), input_stream);
  }

  void
  write_to_zip(const fs::path &output_path, const Output_SDP &sdp,
               const Block_File_Format output_format, const fs::path &temp_dir,
               const size_t num_control_bytes,
               const size_t num_objectives_bytes,
               const std::vector<size_t> &block_info_sizes,
               const std::vector<size_t> &block_data_sizes, Timers &timers)
  {
    // write all files to sdp.zip archive
    Scoped_Timer zip_timer(timers, "zip");
    Archive_Writer writer(output_path);
    // control.json
    {
      Scoped_Timer control_timer(timers, "control");
      auto control_path = get_control_path(temp_dir);
      archive_gzipped_file(control_path, num_control_bytes, writer);
      fs::remove(control_path);
    }
    // objectives.json
    {
      Scoped_Timer objectives_timer(timers, "objectives");
      auto objectives_path = get_objectives_path(temp_dir);
      archive_gzipped_file(objectives_path, num_objectives_bytes, writer);
      fs::remove(objectives_path);
    }

    // block_info_XXX.json
    {
      Scoped_Timer block_info_timer(timers, "block_info");

      // We add block_info before all block_data,
      // so that SDPB can read them fast and skip the rest of archive
      for(size_t block_index = 0; block_index != sdp.num_blocks; ++block_index)
        {
          Scoped_Timer index_timer(timers, std::to_string(block_index));
          const auto block_info_path
            = get_block_info_path(temp_dir, block_index);
          archive_gzipped_file(block_info_path,
                               block_info_sizes.at(block_index), writer);
          fs::remove(block_info_path);
        }
    }
    // block_data_XXX.bin (or .json)
    {
      Scoped_Timer block_data_timer(timers, "block_data");
      for(size_t block_index = 0; block_index != sdp.num_blocks; ++block_index)
        {
          Scoped_Timer index_timer(timers, std::to_string(block_index));
          const auto block_data_path
            = get_block_data_path(temp_dir, block_index, output_format);
          archive_gzipped_file(block_data_path,
                               block_data_sizes.at(block_index), writer);
          fs::remove(block_data_path);
        }
    }
  }

  void check_file_size(const fs::path &file_path, const size_t expected_size)
  {
    const auto size = file_size(file_path);
    ASSERT(size == expected_size, "File size is ", size, ", expected ",
           expected_size, ": ", file_path);
  }

  void check_sdp_directory(const fs::path &temp_dir, size_t num_blocks,
                           Block_File_Format output_format,
                           const size_t num_control_bytes,
                           const size_t num_objectives_bytes,
                           const std::vector<size_t> &block_info_sizes,
                           const std::vector<size_t> &block_data_sizes,
                           Timers &timers)
  {
    Scoped_Timer timer(timers, "check_sdp_dir");
    if(El::mpi::Rank() != 0)
      return;

    // Check
    {
      Scoped_Timer file_count_timer(timers, "file_count");
      size_t file_count = 0;
      for([[maybe_unused]] const auto &it : fs::directory_iterator(temp_dir))
        {
          ++file_count;
        }
      // control.json + objectives.json + (block_info_XXX + block_data_XXX)
      ASSERT(file_count == 2 + 2 * num_blocks, temp_dir, " contains ",
             file_count, " files, ", 2 + 2 * num_blocks, " expected");
    }

    {
      Scoped_Timer file_sizes_timer(timers, "file_sizes");
      check_file_size(get_control_path(temp_dir), num_control_bytes);
      check_file_size(get_objectives_path(temp_dir), num_objectives_bytes);
      for(size_t i = 0; i < num_blocks; ++i)
        {
          check_file_size(get_block_info_path(temp_dir, i),
                          block_info_sizes.at(i));
          check_file_size(get_block_data_path(temp_dir, i, output_format),
                          block_data_sizes.at(i));
        }
    }
  }
}

void write_sdp(const fs::path &output_path, const Output_SDP &sdp,
               Block_File_Format block_file_format, bool zip, Timers &timers,
               const bool debug)
{
  Scoped_Timer write_timer(timers, "write_sdp");

  fs::path temp_dir(output_path);
  temp_dir += "_temp";
  fs::create_directories(temp_dir);
  // We use size_t rather than std::streamsize because MPI treats
  // std::streamsize as an MPI_LONG_INT and then can not MPI_Reduce
  // over it.
  std::vector<size_t> block_info_sizes(sdp.num_blocks, 0);
  std::vector<size_t> block_data_sizes(sdp.num_blocks, 0);

  {
    Scoped_Timer block_files_timer(timers, "block_files");
    for(auto &group : sdp.dual_constraint_groups)
      {
        size_t width = group.constraint_matrix.Width();
        size_t dual_objective_size = sdp.dual_objective_b.size();
        ASSERT(width == dual_objective_size, "width=", width,
               " dual_objective_size=", dual_objective_size);
        {
          Scoped_Timer block_info_timer(
            timers, "block_info_" + std::to_string(group.block_index));
          const auto block_info_path
            = get_block_info_path(temp_dir, group.block_index);
          block_info_sizes.at(group.block_index) = write_data_and_count_bytes(
            block_info_path,
            [&](std::ostream &os) { write_block_info_json(os, group); }, zip);
        }

        {
          Scoped_Timer block_data_timer(
            timers, "block_data_" + std::to_string(group.block_index));
          const auto block_data_path = get_block_data_path(
            temp_dir, group.block_index, block_file_format);
          block_data_sizes.at(group.block_index) = write_data_and_count_bytes(
            block_data_path,
            [&](std::ostream &os) {
              write_block_data(os, group, block_file_format);
            },
            zip, block_file_format == bin);
        }
      }
  }
  {
    Scoped_Timer reduce_timer(timers, "mpi_reduce_block_sizes");

    // for validation:
    auto block_info_sizes_max = block_info_sizes;

    El::mpi::Reduce(block_info_sizes.data(), block_info_sizes.size(),
                    El::mpi::SUM, 0, El::mpi::COMM_WORLD);
    El::mpi::Reduce(block_data_sizes.data(), block_data_sizes.size(),
                    El::mpi::SUM, 0, El::mpi::COMM_WORLD);

    // Validate sizes.
    // For each block, block_info_sizes[block] is positive on exaclty one rank
    // and equals zero on all other ranks.
    // Thus, MPI_Reduce block_info_sizes with SUM and with MAX
    // should yield the same result.
    El::mpi::Reduce(block_info_sizes_max.data(), block_info_sizes_max.size(),
                    El::mpi::MAX, 0, El::mpi::COMM_WORLD);
    if(El::mpi::Rank() == 0)
      {
        for(size_t block_index = 0; block_index < block_info_sizes.size();
            ++block_index)
          {
            ASSERT(block_info_sizes.at(block_index)
                     == block_info_sizes_max.at(block_index),
                   "block_", block_index,
                   " has been processed by more than one rank!");
            ASSERT(block_info_sizes.at(block_index) != 0, "block_info_",
                   block_index, " size is zero");
            ASSERT(block_data_sizes.at(block_index) != 0, "block_data_",
                   block_index, " size is zero");
          }
      }
  }
  const int rank = El::mpi::Rank();
  if(debug)
    {
      print_matrix_sizes(rank, sdp.dual_objective_b,
                         sdp.dual_constraint_groups);
    }
  if(rank == 0)
    {
      // write control.json and objectives.json
      size_t num_control_bytes = write_data_and_count_bytes(
        get_control_path(temp_dir),
        [&](std::ostream &os) {
          write_control_json(os, sdp.num_blocks, sdp.command_arguments);
        },
        zip);
      size_t num_objectives_bytes = write_data_and_count_bytes(
        get_objectives_path(temp_dir),
        [&](std::ostream &os) {
          write_objectives_json(os, sdp.objective_const, sdp.dual_objective_b);
        },
        zip);

      if(zip)
        {
          write_to_zip(output_path, sdp, block_file_format, temp_dir,
                       num_control_bytes, num_objectives_bytes,
                       block_info_sizes, block_data_sizes, timers);
          // Do not call remove_all() to ensure that we
          // don't remove anything useful.
          // This function will fail if temp_dir is not empty.
          fs::remove(temp_dir);
        }
      else
        {
          check_sdp_directory(temp_dir, sdp.num_blocks, block_file_format,
                              num_control_bytes, num_objectives_bytes,
                              block_info_sizes, block_data_sizes, timers);
          fs::rename(temp_dir, output_path);
        }
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
      size_t num_points = group.num_points; // see write_blocks()

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
  ASSERT(my_P * N == my_constraint_matrix_elements,
         "sum(P'*N) != #(B bands): P'=", my_P, ", N=", N,
         "#(B bands)=", my_constraint_matrix_elements);

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
      El::Output("NB: Q is copied over each core group "
                 "(run SDPB with --verbosity=2 to see "
                 "block distribution among core groups).");
      El::Output("---------------------");
    }
}
