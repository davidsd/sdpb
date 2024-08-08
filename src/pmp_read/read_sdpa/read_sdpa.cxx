#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "pmp_read/PMP_File_Parse_Result.hxx"

#include <El.hpp>
#include <cstdlib>
#include <fstream>

namespace
{
  constexpr auto max_size = std::numeric_limits<std::streamsize>::max();

  // https://github.com/vsdp/SDPLIB

  // 1. Comments. The file can begin with arbitrarily many lines of comments. Each line of comments must begin with " or *.
  void read_comments(std::istream &is)
  {
    int next = is.peek();
    while(next == '"' || next == '*')
      {
        is.ignore(max_size, '\n');
        next = is.peek();
        ASSERT(is.good());
      }
  }

  // 2. The first line after the comments contains m, the number of constraint matrices.
  // Additional text on this line after m is ignored.
  void read_m(std::istream &is, size_t &m)
  {
    is >> m;
    is.ignore(max_size, '\n');
    ASSERT(is.good());
  }

  // 3. The second line after the comments contains nblocks, the number of blocks in the block diagonal structure of the matrices.
  // Additional text on this line after nblocks is ignored.
  void read_nblocks(std::istream &is, size_t &nblocks)
  {
    is >> nblocks;
    ASSERT(is.good());
    is.ignore(max_size, '\n');
  }

  // 4. The third line after the comments contains a vector of numbers that give the sizes of the individual blocks.
  // The special characters ',', '(', ')', '{', and '}' can be used as punctuation and are ignored.
  // Negative numbers may be used to indicate that a block is actually a diagonal submatrix.
  // Thus a block size of -5 indicates a 5 by 5 block in which only the diagonal elements are nonzero.
  // TODO: split a diagonal block into 1x1 blocks?
  void read_block_sizes(std::istream &is, std::vector<size_t> &block_sizes,
                        const size_t &nblocks)
  {
    block_sizes.clear();
    block_sizes.reserve(nblocks);
    while(true)
      {
        const int next = is.peek();
        if(next == '\n')
          {
            is.ignore(1);
            return;
          }
        if(next == ',' || next == '(' || next == ')' || next == '{'
           || next == '}' || std::iswspace(next))
          {
            is.ignore(1);
            continue;
          }
        int size;
        is >> size;
        block_sizes.push_back(std::abs(size));
        ASSERT(is.good());
        // In SDPA example.dat-s, there is a comment at the end of comments line, which should be ignored.
        if(block_sizes.size() == nblocks)
          {
            is.ignore(max_size, '\n');
            break;
          }
      }
    ASSERT(is.good());
  }

  // 5. The fourth line after the comments contains the objective function vector c.
  void read_c(std::istream &is, std::vector<El::BigFloat> &c)
  {
    while(true)
      {
        const int next = is.peek();
        if(next == '\n')
          {
            is.ignore(1);
            return;
          }
        if(next == ',' || next == '(' || next == '{' || next == '}'
           || std::iswspace(next))
          {
            is.ignore(1);
            continue;
          }
        El::BigFloat value;
        is >> value;
        c.push_back(value);
      }
    ASSERT(is.good());
  }

  struct SDPA_Matrix_Entry
  {
    // 0..m
    size_t matrix_index{};
    // 1..n_blocks
    size_t block_index{};
    // indices inside block
    size_t i{};
    size_t j{};
    El::BigFloat value;
  };

  // 6. The remaining lines of the file contain entries in the constraint matrices, with one entry per line.
  // The format for each line is
  // <matno> <blkno> <i> <j> <entry>
  // Here <matno> is the number of the matrix to which this entry belongs,
  // <blkno> specifies the block within this matrix,
  // <i> and <j> specify a location within the block,
  // and <entry> gives the value of the entry in the matrix.
  // Note that since all matrices are assumed to be symmetric, only entries in the upper triangle of a matrix are given.
  void read_matrix_entry(std::istream &is, SDPA_Matrix_Entry &entry)
  {
    is >> entry.matrix_index;
    is >> entry.block_index;
    is >> entry.i;
    is >> entry.j;
    is >> entry.value;
    std::ws(is);
    ASSERT(is.eof() || is.good());
  }
}

// Read .dat-s input used in SDPA
// See format description e.g. in https://github.com/vsdp/SDPLIB
PMP_File_Parse_Result
read_sdpa(const std::filesystem::path &input_file,
          const std::function<bool(size_t matrix_index)> &should_parse_matrix)
{
  PMP_File_Parse_Result result;

  std::ifstream is(input_file);
  ASSERT(is.good(), "Could not open ", input_file);

  read_comments(is);

  // m (from SDPA formulation) is equal to N (in SDPB formulation 2.2 or 3.1)
  size_t m;
  read_m(is, m);

  // Number of SDP blocks
  size_t num_blocks;
  read_nblocks(is, num_blocks);
  result.num_matrices = num_blocks;

  std::vector<size_t> block_sizes;
  read_block_sizes(is, block_sizes, num_blocks);
  ASSERT_EQUAL(num_blocks, block_sizes.size());

  std::vector<size_t> block_offsets{0};
  block_offsets.resize(num_blocks);
  for(size_t index = 1; index < num_blocks; ++index)
    {
      block_offsets.at(index)
        = block_offsets.at(index - 1) + block_sizes.at(index - 1);
    }

  std::vector<El::BigFloat> c_objective;
  read_c(is, c_objective);
  ASSERT_EQUAL(m, c_objective.size());

  // Set a_0 = 0, it doesn't affect anything except objective value
  result.objective = {El::BigFloat(0)};
  result.objective->reserve(c_objective.size() + 1);
  for(const auto &value : c_objective)
    {
      // set (a_1..a_N) = -(c_1..c_N)
      // Need to change sign because SDPB formulation requires to maximize a*z,
      // whereas SDPA formulation requires to minimize c*x
      result.objective->push_back(-value);
    }

  // Each matrix element is a vector of degree-0 polynomials: [-F0(i,j), F1(i,j), F2(i,j),...,Fm(i,j)]
  std::vector<std::optional<El::Matrix<Polynomial_Vector>>> pvm_matrices(
    num_blocks);
  for(size_t index = 0; index < pvm_matrices.size(); ++index)
    {
      if(!should_parse_matrix(index))
        continue;
      pvm_matrices.at(index) = El::Matrix<Polynomial_Vector>(
        block_sizes.at(index), block_sizes.at(index));
      // Zero vector of length (m+1)
      Polynomial_Vector zero(m + 1);
      zero.resize(m + 1);
      El::Fill(pvm_matrices.at(index).value(), zero);
    }

  // Fill matrices
  while(!is.eof())
    {
      SDPA_Matrix_Entry entry;
      read_matrix_entry(is, entry);

      if(!should_parse_matrix(entry.block_index - 1))
        continue;

      // entry.i and j and 1-based indices, thus we subtract 1
      int i = entry.i - 1;
      int j = entry.j - 1;
      auto &pvm = pvm_matrices.at(entry.block_index - 1).value();

      auto value = entry.value;
      // The matrix F_0 enters with minus sign in SDPA definition
      if(entry.matrix_index == 0)
        value = -value;

      ASSERT_EQUAL(pvm(i, j).size(), m + 1, DEBUG_STRING(i), DEBUG_STRING(j),
                   DEBUG_STRING(entry.matrix_index),
                   DEBUG_STRING(entry.block_index));
      pvm(i, j).at(entry.matrix_index).coefficients = {value};
      if(i != j)
        pvm(j, i).at(entry.matrix_index).coefficients = {value};
    }

  // Convert matrices to Polynomial_Vector_Matrix
  for(size_t index = 0; index < pvm_matrices.size(); ++index)
    {
      if(!should_parse_matrix(index))
        continue;
      result.parsed_matrices.emplace(
        index, Polynomial_Vector_Matrix(
                 pvm_matrices.at(index).value(), std::nullopt, std::nullopt,
                 std::nullopt, std::nullopt, std::nullopt, std::nullopt));
    }
  return result;
}
