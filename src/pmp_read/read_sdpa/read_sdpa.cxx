#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "pmp_read/PMP_File_Parse_Result.hxx"

#include <El.hpp>
#include <cstdlib>
#include <fstream>
#include <functional>

namespace
{
  enum SDPA_File_Type
  {
    // .dat-s
    sparse,
    // .dat
    dense
  };

  constexpr auto max_size = std::numeric_limits<std::streamsize>::max();

  void skip_line(std::istream &is)
  {
    is.ignore(max_size, '\n');
  }

  bool should_skip(const int token, const bool skip_newline)
  {
    if(!skip_newline && token == '\n')
      return false;
    return token == ',' || token == '(' || token == '{' || token == '}'
           || std::isspace(token);
  }

  void skip_ws(std::istream &is, const bool skip_newline)
  {
    while(should_skip(is.peek(), skip_newline))
      {
        is.ignore(1);
      }
  }

  template <class T> T read_number(std::istream &is, const bool skip_newline)
  {
    T value;
    skip_ws(is, skip_newline);
    is >> value;
    ASSERT(is.good());
    return value;
  }

  // Read vector of a given size from single line.
  // Skip comments at the end of line.
  template <class T>
  std::vector<T> read_vector(std::istream &is, const size_t size)
  {
    std::vector<T> result;
    result.reserve(size);

    for(size_t i = 0; i < size; ++i)
      {
        result.push_back(read_number<T>(is, false));
      }
    skip_line(is);
    ASSERT(is.good());
    return result;
  }

  // https://github.com/vsdp/SDPLIB

  // 1. Comments. The file can begin with arbitrarily many lines of comments. Each line of comments must begin with " or *.
  void read_comments(std::istream &is)
  {
    int next = is.peek();
    while(next == '"' || next == '*')
      {
        skip_line(is);
        next = is.peek();
        ASSERT(is.good());
      }
  }

  // 2. The first line after the comments contains m, the number of constraint matrices.
  // Additional text on this line after m is ignored.
  size_t read_m(std::istream &is)
  {
    size_t m;
    is >> m;
    skip_line(is);
    ASSERT(is.good());
    return m;
  }

  // 3. The second line after the comments contains nblocks, the number of blocks in the block diagonal structure of the matrices.
  // Additional text on this line after nblocks is ignored.
  size_t read_nblocks(std::istream &is)
  {
    size_t nblocks;
    is >> nblocks;
    ASSERT(is.good());
    is.ignore(max_size, '\n');
    return nblocks;
  }

  // 4. The third line after the comments contains a vector of numbers that give the sizes of the individual blocks.
  // The special characters ',', '(', ')', '{', and '}' can be used as punctuation and are ignored.
  // Negative numbers may be used to indicate that a block is actually a diagonal submatrix.
  // Thus a block size of -5 indicates a 5 by 5 block in which only the diagonal elements are nonzero.
  // Further we'll split a diagonal block into 1x1 blocks.
  std::vector<int> read_block_sizes(std::istream &is, const size_t &nblocks)
  {
    return read_vector<int>(is, nblocks);
  }

  // 5. The fourth line after the comments contains the objective function vector c.
  std::vector<El::BigFloat> read_c(std::istream &is, const size_t m)
  {
    return read_vector<El::BigFloat>(is, m);
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
  void read_sparse_matrix_entry(std::istream &is, SDPA_Matrix_Entry &entry)
  {
    is >> entry.matrix_index;
    is >> entry.block_index;
    is >> entry.i;
    is >> entry.j;
    is >> entry.value;
    std::ws(is);
    ASSERT(is.eof() || is.good());
  }

  struct SDPA_Block_Structure
  {
    const size_t dual_dimension;
    const std::vector<int> sdpa_block_sizes;
    std::vector<size_t> sdpb_block_sizes;
    size_t num_sdp_blocks;

  private:
    std::vector<size_t> num_prev_blocks;

  public:
    explicit SDPA_Block_Structure(const size_t m, const size_t nblocks,
                                  const std::vector<int> &block_sizes)
        : dual_dimension(m), sdpa_block_sizes(block_sizes)
    {
      ASSERT_EQUAL(nblocks, block_sizes.size());

      size_t prev_blocks = 0;
      for(const auto &block_size : block_sizes)
        {
          num_prev_blocks.push_back(prev_blocks);
          ASSERT(block_size != 0);
          if(block_size > 0)
            {
              prev_blocks += 1;
              sdpb_block_sizes.push_back(block_size);
            }
          // if block_size is negative, this means a diagonal block.
          // We split it into many 1x1 blocks
          else
            {
              prev_blocks += -block_size;
              sdpb_block_sizes.insert(sdpb_block_sizes.end(), -block_size, 1);
            }
        }
      num_sdp_blocks = prev_blocks;
    }

    void
    get_sdpb_position(const size_t sdpa_block_index, const size_t sdpa_i,
                      const size_t sdpa_j,
                      // output parameters:
                      size_t &sdpb_block_index, size_t &i, size_t &j) const
    {
      const auto &block_size = sdpa_block_sizes.at(sdpa_block_index - 1);
      ASSERT(block_size != 0);
      ASSERT(sdpa_i > 0);
      ASSERT(sdpa_i > 0);
      const auto prev_blocks = num_prev_blocks.at(sdpa_block_index - 1);
      if(block_size > 0)
        {
          sdpb_block_index = prev_blocks;
          i = sdpa_i - 1;
          j = sdpa_j - 1;
        }
      else
        {
          // We represent a diagonal SDP block as a collection of 1x1 blocks.
          ASSERT_EQUAL(sdpa_i, sdpa_j,
                       "Non-diagonal element found in diagonal SDP block!",
                       DEBUG_STRING(sdpa_block_index),
                       DEBUG_STRING(block_size));
          sdpb_block_index = prev_blocks + sdpa_i - 1;
          i = 0;
          j = 0;
        }
    }
  };

  struct SDPB_Matrix_Entry
  {
    // Number of SDP block
    size_t sdp_block_index;
    // Each matrix element is (degree-0) polynomial vector;
    // This is the index inside this vector.
    size_t vector_index;
    size_t i;
    size_t j;
    El::BigFloat value;

    SDPB_Matrix_Entry(const SDPA_Matrix_Entry &entry,
                      const SDPA_Block_Structure &block_structure)
    {
      block_structure.get_sdpb_position(entry.block_index, entry.i, entry.j,
                                        sdp_block_index, i, j);

      this->vector_index = entry.matrix_index;
      // The matrix F_0 enters with minus sign in SDPA definition
      this->value = vector_index == 0 ? -entry.value : entry.value;
    }
  };

  void
  read_sparse_matrices(std::istream &is,
                       const std::function<void(SDPA_Matrix_Entry &&)> &accept)
  {
    // Read all matrix entries in the file
    while(!is.eof())
      {
        SDPA_Matrix_Entry sdpa_entry;
        read_sparse_matrix_entry(is, sdpa_entry);
        std::ws(is);
        accept(std::move(sdpa_entry));
      }
  }

  void
  read_dense_matrices(std::istream &is,
                      const SDPA_Block_Structure &block_structure,
                      const std::function<void(SDPA_Matrix_Entry &&)> &accept)

  {
    // 0..m
    for(size_t matrix_index = 0;
        matrix_index <= block_structure.dual_dimension; ++matrix_index)
      {
        // 1..nblocks
        for(size_t block_index = 1;
            block_index <= block_structure.sdpa_block_sizes.size();
            ++block_index)
          {
            const int sdpa_block_size
              = block_structure.sdpa_block_sizes.at(block_index - 1);
            ASSERT(sdpa_block_size != 0);
            try
              {
                if(sdpa_block_size > 0)
                  {
                    // Parse dense matrix
                    const size_t size = sdpa_block_size;
                    for(size_t i = 1; i <= size; ++i)
                      for(size_t j = 1; j <= size; ++j)
                        {
                          SDPA_Matrix_Entry sdpa_entry;
                          sdpa_entry.block_index = block_index;
                          sdpa_entry.matrix_index = matrix_index;
                          sdpa_entry.i = i;
                          sdpa_entry.j = j;
                          sdpa_entry.value
                            = read_number<El::BigFloat>(is, true);
                          accept(std::move(sdpa_entry));
                        }
                  }
                else
                  {
                    // Parse diagonal matrix entries
                    const size_t size = -sdpa_block_size;
                    for(size_t i = 1; i <= size; ++i)
                      {
                        SDPA_Matrix_Entry sdpa_entry;
                        sdpa_entry.block_index = block_index;
                        sdpa_entry.matrix_index = matrix_index;
                        sdpa_entry.i = i;
                        sdpa_entry.j = i;
                        sdpa_entry.value = read_number<El::BigFloat>(is, true);
                        accept(std::move(sdpa_entry));
                      }
                  }
              }
            catch(const std::exception &e)
              {
                RUNTIME_ERROR("Error when parsing dense SDPA matrix: ",
                              DEBUG_STRING(matrix_index),
                              DEBUG_STRING(block_index),
                              DEBUG_STRING(sdpa_block_size), ":\n", e.what());
              }
          }
      }
  }
}

// Read .dat-s input used in SDPA
// See format description e.g. in https://github.com/vsdp/SDPLIB
PMP_File_Parse_Result
read_sdpa(const std::filesystem::path &input_file,
          const std::function<bool(size_t matrix_index)> &should_parse_matrix)
{
  PMP_File_Parse_Result result;

  const SDPA_File_Type file_type = [&input_file] {
    if(input_file.extension() == ".dat-s")
      return SDPA_File_Type::sparse;
    if(input_file.extension() == ".dat")
      return SDPA_File_Type::dense;
    RUNTIME_ERROR("Unknown SDPA file type, expected .dat or -dat-s: ",
                  input_file.extension());
  }();

  std::ifstream is(input_file);
  ASSERT(is.good(), "Could not open ", input_file);

  read_comments(is);

  const SDPA_Block_Structure block_structure = [&is] {
    // m (from SDPA formulation) is equal to N (in SDPB formulation 2.2 or 3.1)
    const size_t m = read_m(is);

    // Number of SDP blocks
    const size_t num_blocks = read_nblocks(is);

    const auto block_sizes = read_block_sizes(is, num_blocks);
    ASSERT_EQUAL(num_blocks, block_sizes.size());
    return SDPA_Block_Structure(m, num_blocks, block_sizes);
  }();

  result.num_matrices = block_structure.num_sdp_blocks;

  const auto c_objective = read_c(is, block_structure.dual_dimension);
  ASSERT_EQUAL(block_structure.dual_dimension, c_objective.size());

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
  std::vector<std::optional<Simple_Matrix<Polynomial_Vector>>> pvm_matrices(
    block_structure.num_sdp_blocks);
  for(size_t index = 0; index < pvm_matrices.size(); ++index)
    {
      if(!should_parse_matrix(index))
        continue;
      const size_t height = block_structure.sdpb_block_sizes.at(index);
      const size_t width = block_structure.sdpb_block_sizes.at(index);
      pvm_matrices.at(index).emplace(height, width);
      // Zero vector of length (m+1)
      Polynomial_Vector zero(block_structure.dual_dimension + 1);
      for(size_t i = 0; i < height; ++i)
        for(size_t j = 0; j < width; ++j)
          pvm_matrices.at(index).value()(i, j) = zero;
    }

  const auto accept_sdpa_entry = [&](SDPA_Matrix_Entry &&sdpa_entry) {
    const SDPB_Matrix_Entry entry(sdpa_entry, block_structure);
    const auto i = entry.i;
    const auto j = entry.j;

    if(!should_parse_matrix(entry.sdp_block_index))
      return;

    auto &pvm = pvm_matrices.at(entry.sdp_block_index).value();

    ASSERT_EQUAL(pvm(i, j).size(), block_structure.dual_dimension + 1,
                 DEBUG_STRING(i), DEBUG_STRING(j),
                 DEBUG_STRING(sdpa_entry.matrix_index),
                 DEBUG_STRING(sdpa_entry.block_index));
    pvm(i, j).at(entry.vector_index).coefficients = {entry.value};
    if(i != j && file_type == SDPA_File_Type::sparse)
      pvm(j, i).at(entry.vector_index).coefficients = {entry.value};
  };

  switch(file_type)
    {
    case sparse: read_sparse_matrices(is, accept_sdpa_entry); break;
    case dense:
      read_dense_matrices(is, block_structure, accept_sdpa_entry);
      break;
    default: LOGIC_ERROR("Unknown SDPA_File_Type=", file_type);
    }

  skip_ws(is, true);
  if(!is.eof())
    {
      std::string line;
      const auto pos = is.tellg();
      std::getline(is, line);
      RUNTIME_ERROR("Unexpected data at the end of file, position: ", pos,
                    ", next line: ", line);
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
