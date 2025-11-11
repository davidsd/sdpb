#include "Bigint_Syrk_Config.hxx"

// Helper functions. TODO move elsewhere and reuse?
namespace
{
  size_t div_ceil(const size_t a, const size_t b)
  {
    return a / b + (a % b == 0 ? 0 : 1);
  }
  std::vector<size_t> get_input_window_group_heights(
    const std::vector<size_t> &input_height_per_group,
    const size_t split_factor)
  {
    auto heights = input_height_per_group;
    for(auto &height : heights)
      {
        height = div_ceil(height, split_factor);
      }
    return heights;
  }
}

Bigint_Syrk_Config::Bigint_Syrk_Config(
  const El::mpi::Comm &shared_memory_comm, const size_t precision,
  const size_t num_nodes, const std::vector<size_t> &input_height_per_group,
  const size_t input_width, const size_t input_split_factor,
  const size_t output_split_factor)
    : shared_memory_comm(shared_memory_comm),
      precision(precision),
      num_nodes(num_nodes),
      input_height_per_group(input_height_per_group),
      input_height(std::accumulate(input_height_per_group.begin(),
                                   input_height_per_group.end(), 0)),
      input_width(input_width),
      output_dim(input_width),
      input_split_factor(input_split_factor),
      output_split_factor(output_split_factor),
      input_window_height_per_group_per_prime(get_input_window_group_heights(
        input_height_per_group, input_split_factor))
{}
Fmpz_Comb Bigint_Syrk_Config::comb() const
{
  const mp_limb_t bits = precision;
  // Max number of additions during syrk computation
  const size_t max_additions = input_height;
  // Allow negative values
  const int sign = 1;
  return Fmpz_Comb(bits, bits, sign, max_additions);
}
size_t Bigint_Syrk_Config::num_primes() const
{
  return comb().num_primes;
}
size_t Bigint_Syrk_Config::input_window_height() const
{
  return std::accumulate(input_window_height_per_group_per_prime.begin(),
                         input_window_height_per_group_per_prime.end(), 0);
}
size_t Bigint_Syrk_Config::window_width() const
{
  return div_ceil(output_dim, output_split_factor);
}
size_t Bigint_Syrk_Config::input_window_size() const
{
  return input_window_height() * window_width() * num_primes();
}
size_t Bigint_Syrk_Config::num_input_windows() const
{
  return output_split_factor == 1 ? 1 : 2;
}
size_t Bigint_Syrk_Config::output_window_size() const
{
  const size_t dim = window_width();
  return dim * dim * num_primes();
}
size_t Bigint_Syrk_Config::get_reduce_scatter_buffer_size() const
{
  // no reduce-scatter => no RAM needed
  if(num_nodes == 1)
    return 0;

  const size_t dim = div_ceil(output_dim, output_split_factor);

  // Total number of elements in the output submatrix
  // If split factor == 1, we need only half of the matrix.
  // Otherwise, we need the whole matrix for off-diagonal blocks
  const size_t num_elements
    = output_split_factor == 1 ? dim * (dim + 1) / 2 : dim * dim;

  // Each rank needs ~(num_elements/num_ranks) buffer for MPI_Send
  // and the same for MPI_Recv
  return 2 * div_ceil(num_elements, num_nodes);
}
size_t Bigint_Syrk_Config::input_windows_bytes() const
{
  return sizeof(double) * num_input_windows() * input_window_size();
}
size_t Bigint_Syrk_Config::output_window_bytes() const
{
  return sizeof(double) * output_window_size();
}
size_t Bigint_Syrk_Config::node_P_normalizer_bytes() const
{
  return bigfloat_bytes(precision) * input_width * shared_memory_comm.Size();
}
size_t Bigint_Syrk_Config::node_normalize_and_shift_P_bytes() const
{
  // Temporary memory required to AllReduce column norms
  if(shared_memory_comm.Size() == 1)
    return 0;
  // send buffer + recv buffer + extra MPI stuff => factor of 3.
  return 3 * node_P_normalizer_bytes();
}
size_t Bigint_Syrk_Config::node_reduce_scatter_bytes() const
{
  // get_reduce_scatter_buffer_size() returns send_buf + recv_buf.
  // We multiply it by 3/2 to account for extra overhead inside SendRecv.
  return div_ceil(
    bigfloat_bytes(precision) * get_reduce_scatter_buffer_size() * 3, 2);
}
size_t Bigint_Syrk_Config::node_local_bytes() const
{
  const size_t normalizer_size = input_width * shared_memory_comm.Size();
  // Conservative estimate for mpi::AllReduce(column_norms) overhead
  const size_t normalizer_all_reduce_size = 3 * normalizer_size;
  return bigfloat_bytes(precision)
         * (normalizer_size
            + std::max(normalizer_all_reduce_size,
                       get_reduce_scatter_buffer_size()));
}
size_t Bigint_Syrk_Config::node_shmem_bytes() const
{
  // Input and output residue windows
  return input_windows_bytes() + output_window_bytes();
}
size_t Bigint_Syrk_Config::node_total_bytes() const
{
  return node_local_bytes() + node_shmem_bytes();
}
