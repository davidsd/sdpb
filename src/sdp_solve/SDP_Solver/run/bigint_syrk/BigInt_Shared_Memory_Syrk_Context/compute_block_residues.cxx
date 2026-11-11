#include "../BigInt_Shared_Memory_Syrk_Context.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/Fmpz_BigInt.hxx"
#include "sdpb_util/bigint_shared_memory/fmpz/residues.hxx"
#include "sdpb_util/assert.hxx"

#include <optional>

// compute residues and put them to shared window
// NB: input blocks should be BigInt matrix (normalized matrix, multiplied by 2^N)

namespace
{
  // --batchResidues: the residues of a
  // column's entries by Batch_Crt::residues (BLAS, the same doubles as
  // bigint_to_residues) instead of one FLINT multi-mod per entry.
  // Integers |x| < 2^(precision + 2) (normalized and shifted P).
  long batch_residue_bits()
  {
    return long(El::gmp::Precision()) + 2;
  }

  void compute_column_residues_elementwise(
    const El::DistMatrix<El::BigFloat> &block, size_t group_index,
    El::Int residue_row_begin, El::Int global_col, Fmpz_Comb &comb,
    Vertical_Block_Matrix_Residues_Window<double> &block_residues_window,
    Batch_Crt *batch)
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

    // Submatrix is single-column, thus index j is always 0.
    const int j = 0;
    const int jLoc = block_column_submatrix.LocalCol(j);
    if(batch != nullptr)
      {
        const int n = block_column_submatrix.LocalHeight();
        std::vector<const __mpf_struct *> x(n);
        for(int iLoc = 0; iLoc < n; ++iLoc)
          x[iLoc] = block_column_submatrix.GetLocalCRef(iLoc, jLoc)
                      .gmp_float.get_mpf_t();
        std::vector<double> r(size_t(n) * comb.num_primes);
        batch->residues(x.data(), n, batch_residue_bits(), r.data(), 1, n);
        // Write prime by prime: the local rows of one column modulo one
        // prime lie within a few pages of the window, while consecutive
        // primes are prime_stride apart (one page each). Same values as
        // the entry-by-entry order below.
        std::vector<double *> dest(n);
        for(int iLoc = 0; iLoc < n; ++iLoc)
          dest[iLoc] = first_residue_column_submatrix.Buffer(
            block_column_submatrix.GlobalRow(iLoc), j);
        for(size_t k = 0; k < comb.num_primes; ++k)
          {
            const size_t offset = k * block_residues_window.prime_stride;
            const double *rk = r.data() + size_t(n) * k;
            for(int iLoc = 0; iLoc < n; ++iLoc)
              dest[iLoc][offset] = rk[iLoc];
          }
        return;
      }
    Fmpz_BigInt bigint_value;
    for(int iLoc = 0; iLoc < block_column_submatrix.LocalHeight(); ++iLoc)
      {
        const int i = block_column_submatrix.GlobalRow(iLoc);
        bigint_value.from_BigFloat(
          block_column_submatrix.GetLocalCRef(iLoc, jLoc));
        double *data = first_residue_column_submatrix.Buffer(i, j);
        bigint_to_residues(bigint_value, comb, data,
                           block_residues_window.prime_stride);
      }
  }

  // Calculate residues for blocks with consecutive indices
  template <class BlockIterator>
  void compute_column_residues_for_consecutive_blocks(
    El::Int group_index, const BlockIterator &consecutive_blocks_begin,
    const BlockIterator &consecutive_blocks_end, El::Int residue_row_begin,
    El::Int global_col, Fmpz_Comb &comb,
    Vertical_Block_Matrix_Residues_Window<double> &block_residues_window,
    std::vector<double> &column_residues_buffer_temp, Batch_Crt *batch)
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
    if(batch != nullptr)
      {
        std::vector<const __mpf_struct *> x;
        x.reserve(total_height);
        std::for_each(consecutive_blocks_begin, consecutive_blocks_end,
                      [&](const auto &block) {
                        ASSERT(block.IsLocalCol(global_col));
                        ASSERT_EQUAL(block.LocalHeight(), block.Height());
                        const int jLoc = block.LocalCol(global_col);
                        for(int iLoc = 0; iLoc < block.LocalHeight(); ++iLoc)
                          x.push_back(block.GetLocalCRef(iLoc, jLoc)
                                        .gmp_float.get_mpf_t());
                      });
        ASSERT_EQUAL(x.size(), size_t(total_height));
        batch->residues(x.data(), total_height, batch_residue_bits(),
                        column_residues_buffer_temp.data(), 1, prime_stride);
      }
    else
    {
      Fmpz_BigInt bigint_value;
      int data_offset = 0;
      std::for_each(
        consecutive_blocks_begin, consecutive_blocks_end,
        [&](const auto &block) {
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
              bigint_to_residues(bigint_value, comb, data, prime_stride);
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

  using Star_Vr_Block = El::DistMatrix<El::BigFloat, El::STAR, El::VR>;

  void compute_column_residues(
    const size_t group_index,
    const std::vector<El::DistMatrix<El::BigFloat>> &bigint_input_matrix_blocks,
    const El::Int global_col, Fmpz_Comb &comb,
    Vertical_Block_Matrix_Residues_Window<double> &input_block_residues_window,
    std::vector<double> &column_residues_buffer_temp, Batch_Crt *batch,
    const std::vector<std::optional<Star_Vr_Block>> &star_vr_blocks)
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

        // --residuesStarVR: a block on a grid with several rows was
        // copied to [STAR,VR]; the rank owning global_col there writes the
        // whole column at once
        if(!star_vr_blocks.empty() && star_vr_blocks[block_index_local])
          {
            const auto &copy = *star_vr_blocks[block_index_local];
            if(copy.IsLocalCol(global_col))
              compute_column_residues_for_consecutive_blocks(
                group_index, &copy, &copy + 1, residue_row_begin, global_col,
                comb, input_block_residues_window,
                column_residues_buffer_temp, batch);
            residue_row_begin += curr_block->Height();
            continue;
          }

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
            // This will happen if a block is distributed among 4+ ranks.
            // TODO: store blocks in DistMatrix<BigFloat, STAR, VR>
            // instead of default DistMatrix<BigFloat> === DistMatrix<BigFloat, MC, MR>,
            // so that each column will belong to a single rank.
            // See Elemental: A New Framework for Distributed Memory Dense Matrix Computations, Poulsen et al. (2016)
            // https://www.cs.utexas.edu/~flame/pubs/Elemental1.pdf (page 12)
            compute_column_residues_elementwise(
              *curr_block, group_index, residue_row_begin, global_col, comb,
              input_block_residues_window, batch);

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
          column_residues_buffer_temp, batch);

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
  Vertical_Block_Matrix_Residues_Window<double> &grouped_block_residues_window,
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

    // --batchResidues: residues by BLAS (see above)
    Batch_Crt *batch = nullptr;
    {
      if(options.batch_residues && Batch_Crt::supported(comb))
        {
          if(!batch_crt)
            batch_crt = std::make_unique<Batch_Crt>(comb);
          batch = batch_crt.get();
        }
    }
    // --residuesStarVR: a block on a grid
    // with several process rows has no rank holding whole columns, so its
    // residues were written entry by entry, num_primes scattered writes per
    // entry (the dominant cost of compute_and_write on large SDPs). Copy such
    // blocks to [STAR,VR] first (each column on one rank; SDPB's Readme TODO):
    // the same values go to the same places of the window, column by column.
    // The decision uses the grid shape, so it is the same on all ranks of the
    // group (the copy is collective over the block's grid).
    std::vector<std::optional<Star_Vr_Block>> star_vr_blocks;
    {
      if(residues_star_vr)
        {
          Scoped_Timer star_vr_timer(timers, "star_vr");
          star_vr_blocks.resize(block_views.size());
          for(size_t b = 0; b < block_views.size(); ++b)
            if(block_views[b].Grid().Height() > 1)
              star_vr_blocks[b].emplace(block_views[b]);
        }
    }
    std::vector<double> column_residues_buffer_temp;
    for(int global_col = 0; global_col < width; ++global_col)
      {
        compute_column_residues(group_index, block_views, global_col, comb,
                                grouped_block_residues_window,
                                column_residues_buffer_temp, batch,
                                star_vr_blocks);
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