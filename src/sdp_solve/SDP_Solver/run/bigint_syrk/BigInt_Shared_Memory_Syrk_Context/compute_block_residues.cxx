#include "../BigInt_Shared_Memory_Syrk_Context.hxx"
#include "../fmpz/Fmpz_BigInt.hxx"
#include "../fmpz/fmpz_mul_blas_util.hxx"
#include "sdpb_util/assert.hxx"

// compute residues and put them to shared window
// NB: input blocks should be BigInt matrix (normalized matrix, multiplied by 2^N)

namespace
{
  void compute_column_residues_elementwise(
    const El::DistMatrix<El::BigFloat> &block, size_t group_index,
    El::Int residue_row_begin, El::Int global_col, Fmpz_Comb &comb,
    Block_Residue_Matrices_Window<double> &block_residues_window)
  {
    if(!block.IsLocalCol(global_col))
      return;

    // Block submatrix for a given row range and a given column
    const auto block_column_submatrix
      = block(El::Range(0, block.Height()), global_col);
    ASSERT(block_column_submatrix.Viewing());

    // Residues of blocks on a current MPI group
    // modulo the first prime are written here
    auto &group_residues_matrix
      = block_residues_window.block_residues.at(0).at(group_index);

    // Part of shared memory window containing
    // residues of block_column_submatrix modulo the first prime number.
    const El::Range<El::Int> residue_I(residue_row_begin,
                                       residue_row_begin + block.Height());
    const El::Range<El::Int> residue_J(global_col);
    ASSERT(residue_I.beg < group_residues_matrix.Height(),
           DEBUG_STRING(residue_I.beg),
           DEBUG_STRING(group_residues_matrix.Height()));
    ASSERT(residue_I.end <= group_residues_matrix.Height());
    auto first_residue_column_submatrix
      = group_residues_matrix(residue_I, residue_J);
    ASSERT(first_residue_column_submatrix.Viewing());

    Fmpz_BigInt bigint_value;
    // Submatrix is single-column, thus index j is always 0.
    const int j = 0;
    const int jLoc = block_column_submatrix.LocalCol(j);
    for(int iLoc = 0; iLoc < block_column_submatrix.LocalHeight(); ++iLoc)
      {
        const int i = block_column_submatrix.GlobalRow(iLoc);
        bigint_value.from_BigFloat(
          block_column_submatrix.GetLocalCRef(iLoc, jLoc));
        double *data = first_residue_column_submatrix.Buffer(i, j);
        fmpz_multi_mod_uint32_stride(data, block_residues_window.prime_stride,
                                     bigint_value.value, comb);
      }
  }

  // Calculate residues for blocks with consecutive indices
  template <class BlockIterator>
  void compute_column_residues_for_consecutive_blocks(
    El::Int group_index, const BlockIterator &consecutive_blocks_begin,
    const BlockIterator &consecutive_blocks_end, El::Int residue_row_begin,
    El::Int global_col, Fmpz_Comb &comb,
    Block_Residue_Matrices_Window<double> &block_residues_window,
    std::vector<double> &column_residues_buffer_temp)
  {
    auto total_height = std::transform_reduce(
      consecutive_blocks_begin, consecutive_blocks_end, 0, std::plus{},
      [](auto &block) { return block.Height(); });
    // Stack blocks on top of each other,
    // calculate residues for each column locally and then copy them to the window.
    // This should reduce the number of shared memory access operations,
    // which can be ~10x slower than local memory access.

    size_t prime_stride = total_height;

    // (column residues for prime_1),(column residues for prime_2),... stored in a single array.
    column_residues_buffer_temp.resize(total_height * comb.num_primes);

    // Compute locally residues for a given column global_col and all blocks
    {
      Fmpz_BigInt bigint_value;
      int data_offset = 0;
      std::for_each(
        consecutive_blocks_begin, consecutive_blocks_end,
        [&](const El::DistMatrix<El::BigFloat> &block) {
          // For a given column jLoc, calculate all residues locally and put to column_residues
          ASSERT(block.IsLocalCol(global_col));
          ASSERT_EQUAL(block.LocalHeight(), block.Height());
          int jLoc = block.LocalCol(global_col);
          for(int iLoc = 0; iLoc < block.LocalHeight(); ++iLoc)
            {
              bigint_value.from_BigFloat(block.GetLocalCRef(iLoc, jLoc));
              // pointer to the first residue
              double *data
                = column_residues_buffer_temp.data() + data_offset + iLoc;
              fmpz_multi_mod_uint32_stride(data, prime_stride,
                                           bigint_value.value, comb);
            }
          data_offset += block.LocalHeight();
        });
      ASSERT_EQUAL(data_offset, total_height);
    }

    {
      // For each prime, copy column residues to shared memory window
      for(size_t prime_index = 0; prime_index < comb.num_primes; ++prime_index)
        {
          int j = global_col;

          auto &residue_matrix
            = block_residues_window.block_residues.at(prime_index)
                .at(group_index);

          double *src
            = column_residues_buffer_temp.data() + prime_index * prime_stride;
          double *dest = residue_matrix.Buffer(residue_row_begin, j);

          // according to C++ docs, memcpy is the fastest way to copy memory
          std::memcpy(dest, src, sizeof(double) * total_height);
        }
    }
  }

  void compute_column_residues(
    const size_t group_index,
    const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
    const El::Int global_col, Fmpz_Comb &comb,
    Block_Residue_Matrices_Window<double> &input_block_residues_window,
    std::vector<double> &column_residues_buffer_temp)
  {
    // offset for current block in the residues window
    int residue_row_begin = 0;
    // For each column, merge consecutive blocks (e.g. index=3,4,5,6) that own this column locally,
    // compute residues for this column and memcpy the residues array to block_residues_window.
    const size_t num_blocks_local = bigint_input_matrix_blocks.size();
    for(size_t block_index_local = 0; block_index_local < num_blocks_local;
        ++block_index_local)
      {
        const auto curr_block
          = bigint_input_matrix_blocks.begin() + block_index_local;

        if(!curr_block->IsLocalCol(global_col))
          {
            residue_row_begin += curr_block->Height();
            continue;
          }

        if(curr_block->Height() != curr_block->LocalHeight())
          {
            // If the rank owns only part of the column,
            // then residues of local elements are not stored consecutively in input_block_residues_window.
            // In that case, we cannot compute residues for the whole column locally
            // and then memcpy them to the shared window.
            // Thus we have to write residues one by one.
            //
            // This will happen if a block is distributed among 4+ nodes.
            // TODO: store blocks in DistMatrix<BigFloat, STAR, VR>
            // instead of default DistMatrix<BigFloat> === DistMatrix<BigFloat, MC, MR>,
            // so that each column will belong to a single rank.
            // See Elemental: A New Framework for Distributed Memory Dense Matrix Computations, Poulsen et al. (2016)
            // https://www.cs.utexas.edu/~flame/pubs/Elemental1.pdf (page 12)
            compute_column_residues_elementwise(
              *curr_block, group_index, residue_row_begin, global_col, comb,
              input_block_residues_window);

            residue_row_begin += curr_block->Height();
            continue;
          }

        // Take all blocks that:
        // - have consecutive indices
        // - own all elements of the current column global_col
        // Then, for each prime, all residues from these blocks
        // are stored consecutively in input_block_residues_window.
        int total_consecutive_height = curr_block->Height();
        size_t num_consecutive_blocks = 1;
        for(size_t delta = 1; block_index_local + delta < num_blocks_local;
            ++delta)
          {
            const auto next_block = curr_block + delta;
            if(!next_block->IsLocalCol(global_col))
              break;
            if(next_block->LocalHeight() != next_block->Height())
              break;

            ++num_consecutive_blocks;
            total_consecutive_height += next_block->Height();
          }

        compute_column_residues_for_consecutive_blocks(
          group_index, curr_block, curr_block + num_consecutive_blocks,
          residue_row_begin, global_col, comb, input_block_residues_window,
          column_residues_buffer_temp);

        residue_row_begin += total_consecutive_height;

        // Skip blocks processed in the previous line.
        block_index_local += num_consecutive_blocks - 1;
        // block_index_local now points to the last processed block,
        // i.e. at the start of the next iteration it will point to the first unprocessed block.
      }
  }
}

// Each rank computes residues for locally stored block elements and puts them into shared memory window.
// Memory access to the window is more expensive (up to ~10x in our tests) than access to local memory.
// Because of that, we want to compute many residues locally and then memcpy the resulting array to the window.
// For each prime, residues from all blocks on a node are stored in a Matrix<double> in a column-major order.
// It means that, for each column, elements from all blocks are stored consecutively.
//
// If a rank owns, e.g., the whole column 3 for blocks with indices (2,3,4,5),
// Then we can calculate residues for this column from all blocks combined, write them to a local array
// and then memcpy this array to the shared memory window.
// See compute_column_residues_for_consecutive_blocks().
//
// If residues don't form a contiguous memory range (e.g. rank owns only half of the column elements),
// then we calculate residues for each element separately and put them into shared memory window one by one.
// See compute_column_residues_elementwise().
//
// We calculate first skip_rows rows,
// and then fill the input window with as many remaining rows as we can
// (if the window height is small, we have to call compute_block_residues() several times)
void BigInt_Shared_Memory_Syrk_Context::compute_block_residues(
  Block_Residue_Matrices_Window<double> &grouped_block_residues_window,
  const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
  El::Int skip_rows, El::Range<El::Int> col_range, Timers &timers,
  El::Matrix<int32_t> &block_timings_ms)
{
  Scoped_Timer compute_residues_timer(timers, "compute_residues");
  {
    Scoped_Timer compute_and_write_timer(timers, "compute_and_write");
    // How many block rows we want to write to the window
    int target_height = input_group_height_per_prime();
    int width = col_range.end - col_range.beg;
    ASSERT(col_range.beg >= 0 && width > 0, DEBUG_STRING(col_range.beg),
           DEBUG_STRING(col_range.end));

    // block_views contains subset of block rows that will be processed
    // in the current step.
    std::vector<El::DistMatrix<El::BigFloat>> block_views;
    for(const auto &block : bigint_input_matrix_blocks)
      {
        if(target_height == 0)
          break;
        ASSERT(target_height > 0, DEBUG_STRING(target_height));
        const auto height = block.Height();
        ASSERT(skip_rows >= 0, DEBUG_STRING(skip_rows));
        if(skip_rows >= height)
          {
            skip_rows -= height;
            continue;
          }
        const int row_begin = skip_rows;
        const int row_end = std::min(height, row_begin + target_height);

        auto &block_view = block_views.emplace_back();
        El::LockedView(block_view, block, El::Range<int>(row_begin, row_end),
                       col_range);
        ASSERT(block_view.Viewing());

        skip_rows = 0;
        target_height -= row_end - row_begin;
      }
    if(target_height > 0)
      {
        // Fill remaining rows with zeros.
        // TODO: write zeros directly to residue window?
        // That would be faster, but it's hardly a bottleneck.
        // So we go with a simple solution for now.
        const auto &grid = bigint_input_matrix_blocks.at(0).Grid();
        auto &zero_block
          = block_views.emplace_back(target_height, width, grid);
        El::Zero(zero_block);
        target_height = 0;
      }

    {
      // Sanity check
      const int block_views_height = std::transform_reduce(
        block_views.begin(), block_views.end(), 0, std::plus{},
        [](const auto &matrix) { return matrix.Height(); });
      ASSERT_EQUAL(block_views_height, input_group_height_per_prime());
    }

    std::vector<double> column_residues_buffer_temp;
    for(int global_col = 0; global_col < width; ++global_col)
      {
        compute_column_residues(group_index, block_views, global_col, comb,
                                grouped_block_residues_window,
                                column_residues_buffer_temp);
      }

    // Update block timings.
    // Take total time and split it among local blocks, proportionally to block height

    Scoped_Timer block_timing_timer(timers, "block_timing");
    const auto total_time_ms = compute_and_write_timer.elapsed_milliseconds();

    ASSERT_EQUAL(bigint_input_matrix_blocks.size(),
                 block_index_local_to_global.size());
    constexpr auto get_size
      = [](auto &block) { return block.LocalHeight() * block.LocalWidth(); };
    const auto total_size = std::transform_reduce(
      bigint_input_matrix_blocks.cbegin(), bigint_input_matrix_blocks.cend(),
      0, std::plus{}, get_size);
    // Average time spent on one block element
    const double time_per_element = (double)total_time_ms / total_size;

    // total_size=0 e.g. if we have single 1x1 block distributed to 2 ranks.
    // Then rank=1 doesn't have any elements
    if(total_size > 0)
      {
        for(size_t local_block_index = 0;
            local_block_index < block_index_local_to_global.size();
            ++local_block_index)
          {
            auto block_size
              = get_size(bigint_input_matrix_blocks.at(local_block_index));
            double time = block_size * time_per_element;
            auto global_block_index
              = block_index_local_to_global.at(local_block_index);
            block_timings_ms(global_block_index, 0) += std::round(time);
          }
      }
  }

  // wait for all ranks to fill input_block_residues_window
  Scoped_Timer fence_timer(timers, "fence");
  grouped_block_residues_window.Fence();
}