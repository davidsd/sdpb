#include "memory_estimates.hxx"

#include "Memory_Tracker.hxx"
#include "Proc_Meminfo.hxx"
#include "assert.hxx"
#include "ostream/pretty_print_bytes.hxx"

#include <iomanip>
#include <unistd.h>

size_t bigfloat_bytes()
{
  return sizeof(El::BigFloat) + El::gmp::num_limbs * sizeof(mp_limb_t);
}
size_t bigfloat_bytes(const mp_bitcnt_t precision)
{
  // Copied from El::SetPrecision
  const auto num_limbs = (std::max(precision, static_cast<mp_bitcnt_t>(53))
                          + 2 * GMP_NUMB_BITS - 1)
                           / GMP_NUMB_BITS
                         + 1;
  return sizeof(El::BigFloat) + num_limbs * sizeof(mp_limb_t);
}

// TODO rewrite using --maxMemory limit
size_t get_max_shared_memory_bytes(
  const size_t nonshared_memory_required_per_node_bytes,
  const Environment &env, const Verbosity verbosity)
{
  const size_t mem_total_bytes = env.node_mem_total();
  if(mem_total_bytes == 0)
    return 0;

  // Total memory required for all ranks on a node

  size_t max_shared_memory_bytes = 0;
  if(env.comm_shared_mem.Rank() == 0)
    {
      std::ostringstream ss;
      El::BuildStream(ss, "node=", env.node_index(), ": ");

      El::BuildStream(
        ss, "\n\tMemTotal: ", pretty_print_bytes(mem_total_bytes, false),
        "\n\tRequired memory estimate (excluding shared memory windows): ",
        pretty_print_bytes(nonshared_memory_required_per_node_bytes, false));

      if(nonshared_memory_required_per_node_bytes > mem_total_bytes)
        {
          // This is certainly not enough, but at least
          // we'll print sizes in BigInt_Shared_Memory_Syrk_Context
          max_shared_memory_bytes = 0.5 * mem_total_bytes;
          El::BuildStream(
            ss,
            "\n\tSDPB will set --maxSharedMemory to 50% of MemTotal, i.e. ",
            pretty_print_bytes(max_shared_memory_bytes, false));
          El::BuildStream(ss, "\n\tSDPB will probably fail with OOM. Consider "
                              "increasing number of nodes or RAM per node.");
          if(verbosity < Verbosity::debug)
            {
              El::BuildStream(ss, "\n\tTo print detailed memory estimates, "
                                  "run SDPB with --verbosity debug.");
            }
        }
      else
        {
          // ad-hoc coefficient 0.5 to leave some free RAM
          max_shared_memory_bytes
            = 0.5
              * (mem_total_bytes - nonshared_memory_required_per_node_bytes);
          El::BuildStream(ss,
                          "\n\tTo prevent OOM, "
                          "SDPB will set --maxSharedMemory to 50% of the "
                          "remaining memory, i.e. ",
                          pretty_print_bytes(max_shared_memory_bytes, false));
          El::BuildStream(ss,
                          "\n\tIn case of OOM, consider increasing number of "
                          "nodes and/or decreasing --maxSharedMemory limit.");
        }
      if(verbosity >= Verbosity::regular)
        {
          if(env.node_index() == 0 || verbosity >= Verbosity::debug)
            PRINT_WARNING(ss.str());
        }
    }
  // All ranks on a node should have the same limit
  El::mpi::Broadcast(max_shared_memory_bytes, 0, env.comm_shared_mem);

  return max_shared_memory_bytes;
}
size_t get_heap_allocated_bytes(const El::BigFloat &)
{
  return bigfloat_bytes();
}
size_t get_heap_allocated_bytes(const El::AbstractDistMatrix<El::BigFloat> &m)
{
  return m.AllocatedMemory() * bigfloat_bytes();
}
size_t get_heap_allocated_bytes(const El::Matrix<El::BigFloat> &m)
{
  return m.MemorySize() * bigfloat_bytes();
}
void print_allocation_message_per_node(const Environment &env,
                                       const std::string &name, size_t bytes)
{
  bytes = El::mpi::Reduce(bytes, 0, env.comm_shared_mem);
  if(env.comm_shared_mem.Rank() == 0)
    {
      El::Output("node=", env.node_index(), ": allocate ", name, ": ",
                 pretty_print_bytes(bytes));
    }
}
size_t get_trsm_bytes(const int height, const int width, const int grid_height,
                      const int grid_width)
{
  ASSERT(height > 0, DEBUG_STRING(height));
  ASSERT(width > 0, DEBUG_STRING(width));
  ASSERT(grid_height > 0, DEBUG_STRING(grid_height));
  ASSERT(grid_width > 0, DEBUG_STRING(grid_width));

  const int m = height;
  const int n = width;
  const int grid_size = grid_height * grid_width;
  // algorithmic block size in Elemental, default value: 128
  const int bsize = El::Blocksize();
  const size_t L = m * m;
  const size_t X = m * n;

  // Maximum size of temporary matrices
  // across all iterations in LLNLarge() or LLNMedium()
  // One can check that it's always the first iteration, k = 0.
  size_t all_matrices_size = 0;
  // Estimate memory used for MPI communication
  size_t communication_size = 0;

  const int nb = El::Min(bsize, m);
  const int ind1 = nb;
  const int ind2 = m - nb;
  // Views over submatrices of L and X
  const size_t L11 = ind1 * ind1;
  const size_t L21 = ind2 * ind1;
  const size_t X1 = ind1 * n;
  // const size_t X2 = ind2 * n;

  // DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
  // (makes a copy of L)
  const size_t LProx = L;
  // DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
  // (makes a copy of X)
  const size_t XProx = X;
  all_matrices_size += LProx + XProx;

  // Temporary matrices
  // DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
  // Data is duplicated over all grid elements
  // To learn more about different DistMatrix distributions, see
  // https://bootstrapcollaboration.gitlab.io/elemental-web/documentation/dev/core/dist_matrix/DM.html
  const size_t L11_STAR_STAR = L11 * grid_size;
  // DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
  // Data is duplicated over all grid columns
  const size_t L21_MC_STAR = L21 * grid_width;
  all_matrices_size += L11_STAR_STAR + L21_MC_STAR;

  // Choose the largest matrix that require communication
  // for copying from/to (MC,MR) distribution
  if(grid_size > 1)
    communication_size = L11_STAR_STAR;
  if(grid_height > 1)
    communication_size = std::max(communication_size, L21_MC_STAR);

  // LLNLarge
  // https://gitlab.com/bootstrapcollaboration/elemental/-/blob/3ee9f404e632d9d0a4751e9b6717da51e59b5697/src/blas_like/level3/Trsm/LLN.hpp#L19
  if(n > 5 * grid_size)
    {
      // DistMatrix<F,STAR,VR  > X1_STAR_VR(g);
      // Data is not duplicated
      const size_t X1_STAR_VR = X1;
      // DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
      // Data is duplicated over all grid rows
      const size_t X1_STAR_MR = X1 * grid_height;

      all_matrices_size += X1_STAR_VR + X1_STAR_MR;

      if(grid_height > 1)
        {
          communication_size
            = std::max({communication_size, X1_STAR_VR, X1_STAR_MR});
        }
    }
  // LLNMedium
  else
    {
      // DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);
      // Data is duplicated over all grid rows
      const size_t X1Trans_MR_STAR = X1 * grid_height;
      all_matrices_size += X1Trans_MR_STAR;
      if(grid_height > 1)
        {
          communication_size = std::max(communication_size, X1Trans_MR_STAR);
        }
    }

  // Factor of 3 before communication_size is somewhat arbitrary.
  // It seems to work well for our experimental data
  // for wide matrices and large grid sizes.
  // Checked on Expanse HPC e.g. for:
  // height = [31,62,125], width = 125*[427,854,1708], precision=448, grid_size=1..64
  return bigfloat_bytes() * (all_matrices_size + 3 * communication_size);
}
size_t get_trmm_bytes(const int height, const int width, const int grid_height,
                      const int grid_width)
{
  ASSERT(height > 0, DEBUG_STRING(height));
  ASSERT(width > 0, DEBUG_STRING(width));
  ASSERT(grid_height > 0, DEBUG_STRING(grid_height));
  ASSERT(grid_width > 0, DEBUG_STRING(grid_width));

  const int m = height;
  const int n = width;
  const int grid_size = grid_height * grid_width;
  // algorithmic block size in Elemental, default value: 128
  const int bsize = El::Blocksize();
  const size_t X = m * n;
  const size_t L = n * n;

  // Maximum size of temporary matrices
  // across all iterations in RLNA() or RLNC()
  // One can check that it's always the first iteration, k = 0.
  size_t all_matrices_size = 0;
  // Estimate memory used for MPI communication
  size_t communication_size = 0;

  // DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
  // (makes a copy of L)
  const size_t LProx = L;
  // DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
  // (makes a copy of X)
  const size_t XProx = X;
  all_matrices_size += LProx + XProx;

  // RLNA
  // https://gitlab.com/bootstrapcollaboration/elemental/-/blob/3ee9f404e632d9d0a4751e9b6717da51e59b5697/src/blas_like/level3/Trmm/RLN.hpp#L72
  if(width > 5 * height)
    {
      const int nb = El::Min(bsize, m);
      const size_t X1 = nb * n;
      const size_t Z1 = n * nb;

      // DistMatrix<T,STAR,VC  > X1_STAR_VC(g);
      // Data is not duplicated
      const size_t X1_STAR_VC = X1;
      // DistMatrix<T,STAR,MC  > X1_STAR_MC(g);
      // Data is duplicated over all grid columns
      const size_t X1_STAR_MC = X1 * grid_width;
      // DistMatrix<T,MR,  STAR> Z1Trans_MR_STAR(g);
      // Data is duplicated over all grid rows
      const size_t Z1Trans_MR_STAR = Z1 * grid_height;
      // DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);
      // Data is not duplicated
      const size_t Z1Trans_MR_MC = Z1;

      all_matrices_size
        += X1_STAR_VC + X1_STAR_MC + Z1Trans_MR_STAR + Z1Trans_MR_MC;

      if(grid_height > 1)
        {
          communication_size = std::max(communication_size, Z1Trans_MR_STAR);
        }
      if(grid_width > 1)
        {
          communication_size = std::max(communication_size, X1_STAR_MC);
        }
    }
  // RLNC
  // https://gitlab.com/bootstrapcollaboration/elemental/-/blob/3ee9f404e632d9d0a4751e9b6717da51e59b5697/src/blas_like/level3/Trmm/RLN.hpp#L173
  else
    {
      size_t RNLC_matrices_size = 0;
      size_t RNLC_communication_size = 0;
      // In this case, we don't know which iterations has max memory usage.
      // It should be one of the two last iterations, depending on (n % bsize).
      // TODO: check only them?
      // This shouldn't be a bottleneck, so let's keep the explicit for loop.
      for(int k = 0; k < n; k += bsize)
        {
          size_t curr_matrices_size = 0;
          size_t curr_communication_size = 0;
          const int nb = El::Min(bsize, n - k);
          const size_t L10 = nb * k;
          const size_t L11 = nb * nb;
          // const size_t X0 = m * k;
          const size_t X1 = m * nb;

          // DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
          // Data is duplicated over all grid elements
          const size_t L11_STAR_STAR = L11 * grid_height * grid_width;
          // DistMatrix<T,MR,  STAR> L10Trans_MR_STAR(g);
          // Data is duplicated over all grid rows
          const size_t L10Trans_MR_STAR = L10 * grid_height;
          // DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
          // Data is not duplicated
          const size_t X1_VC_STAR = X1;
          // DistMatrix<T,MC,  STAR> X1_MC_STAR(g);
          // Data is duplicated over all grid columns
          const size_t X1_MC_STAR = X1 * grid_width;

          curr_matrices_size
            += L11_STAR_STAR + L10Trans_MR_STAR + X1_VC_STAR + X1_MC_STAR;
          if(grid_height > 1)
            {
              curr_communication_size = std::max(
                {curr_communication_size, L11_STAR_STAR, L10Trans_MR_STAR});
            }
          if(grid_width > 1)
            {
              curr_communication_size
                = std::max(curr_communication_size, X1_MC_STAR);
            }
          RNLC_matrices_size
            = std::max(RNLC_matrices_size, curr_matrices_size);
          RNLC_communication_size
            = std::max(RNLC_communication_size, curr_communication_size);
        }
      all_matrices_size += RNLC_matrices_size;
      communication_size
        = std::max(communication_size, RNLC_communication_size);
    }

  // Factor of 3 before communication_size is somewhat arbitrary.
  // See comments in get_trsm_bytes() above.
  return bigfloat_bytes() * (all_matrices_size + 3 * communication_size);
}

size_t get_cholesky_bytes(const El::UpperOrLower uplo, int height, int width,
                          int grid_height, int grid_width)
{
  using Allocation = Memory_Tracker::Allocation;
  using Scope = Memory_Tracker::Scope;

  const auto grid_size = grid_height * grid_width;

  // Approximate factor for redistributing a matrix of size X
  // for MPI communication, e.g. copying a matrix of size X
  // to another rank should require ~ 3 * X memory
  constexpr int mpi_copy_factor = 3;

  // Reproduce El::Cholesky::LowerVariant3Blocked()
  // https://gitlab.com/bootstrapcollaboration/elemental/-/blob/6c11ea5df08b9f37193d7eb1df2400ed4a3a4d28/src/lapack_like/factor/Cholesky/LowerVariant3.hpp#L76
  Memory_Tracker t("cholesky");
  {
    using El::Int, El::Min, El::Range;
    ASSERT_EQUAL(height, width);
    Allocation AProx(t, "AProx", height * width * bigfloat_bytes());
    const Int n = height;
    const Int bsize = El::Blocksize();
    const Int STAR_STAR = grid_size;
    const Int VC_STAR = 1;
    const Int VR_STAR = 1;
    const Int MC_STAR = grid_width;
    const Int MR_STAR = grid_height;
    const Int STAR_VC = 1;
    const Int STAR_VR = 1;
    const Int STAR_MC = grid_width;
    const Int STAR_MR = grid_height;

    for(Int k = 0; k < n; k += bsize)
      {
        const Int nb = Min(bsize, n - k);
        const Int ind1 = nb;
        const Int ind2 = n - (k + nb);
        const auto A11_bytes = ind1 * ind1 * bigfloat_bytes();
        const auto A12_bytes = ind1 * ind2 * bigfloat_bytes();
        const auto A21_bytes = ind2 * ind1 * bigfloat_bytes();
        const auto A22_bytes = ind2 * ind2 * bigfloat_bytes();

        Allocation A11_STAR_STAR(t, "A11_STAR_STAR", A11_bytes * STAR_STAR);

        if(uplo == El::UPPER)
          {
            Allocation A12_STAR_VR(t, "A12_STAR_VR", A21_bytes * STAR_VR);
            Allocation A12_STAR_MC(t, "A12_STAR_MC", A21_bytes * STAR_MC);
            Allocation A12_STAR_MR(t, "A12_STAR_MR", A21_bytes * STAR_MR);
            if(grid_size > 1)
              {
                // Memory footprint for different operations involving MPI communication.
                std::ignore
                  = Allocation(t, "A11_STAR_STAR = A11;",
                               mpi_copy_factor * A11_bytes * grid_size);
                std::ignore = Allocation(t, "A12_STAR_VR = A12;",
                                         mpi_copy_factor * A12_bytes);
                std::ignore
                  = Allocation(t, "A12_STAR_MC = A12_STAR_VR;",
                               mpi_copy_factor * A12_bytes * STAR_MC);
                std::ignore
                  = Allocation(t, "A12_STAR_MR = A12_STAR_VR;",
                               mpi_copy_factor * A12_bytes * STAR_MR);
                std::ignore
                  = Allocation(t, "A12 = A12_STAR_MR;",
                               mpi_copy_factor * A12_bytes * STAR_MR);
              }
          }
        else if(uplo == El::LOWER)
          {
            Allocation A21_VC_STAR(t, "A21_VC_STAR", A21_bytes * VC_STAR);

            Allocation A21_VR_STAR(t, "A21_VR_STAR", A21_bytes * VR_STAR);

            Allocation A21Trans_STAR_MC(t, "A21Trans_STAR_MC",
                                        A21_bytes * STAR_MC);
            Allocation A21Adj_STAR_MR(t, "A21Adj_STAR_MR",
                                      A21_bytes * STAR_MR);
            if(grid_size > 1)
              {
                // Memory footprint for different operations involving MPI communication.
                std::ignore
                  = Allocation(t, "A11_STAR_STAR = A11;",
                               mpi_copy_factor * A11_bytes * grid_size);
                std::ignore = Allocation(t, "A21_VC_STAR = A21;",
                                         mpi_copy_factor * A21_bytes);
                std::ignore = Allocation(t, "A21_VR_STAR = A21_VC_STAR;",
                                         mpi_copy_factor * A21_bytes);
                {
                  Scope Transpose_A21_VC_STAR(
                    t, "Transpose( A21_VC_STAR, A21Trans_STAR_MC );");
                  // see El::transpose::PartialColAllGather
                  Allocation ATrans(t, "ATrans", A21_bytes * grid_width);
                  std::ignore = Allocation(
                    t, "Copy()", mpi_copy_factor * A21_bytes * grid_width);
                }
                {
                  Scope Adjoint_A21_VR_STAR(
                    t, "Adjoint( A21_VR_STAR, A21Adj_STAR_MR );");
                  // see El::transpose::PartialColAllGather
                  Allocation ATrans(t, "ATrans", A21_bytes * grid_height);
                  std::ignore = Allocation(
                    t, "Copy()", mpi_copy_factor * A21_bytes * grid_height);
                }
                // El::transpose::RowFilter creates one copy of A21Trans_STAR_MC
                std::ignore
                  = Allocation(t, "Transpose( A21Trans_STAR_MC, A21 );",
                               A21_bytes * grid_width);
              }
          }
        else
          {
            LOGIC_ERROR(DEBUG_STRING(uplo));
          }
      }
  }
  return t.peak_memory();
}
size_t shmem_overhead_bytes(const El::mpi::Comm &comm_shared_mem,
                            const size_t shmem_bytes)
{
  const size_t num_ranks = comm_shared_mem.Size();

  long int page_size = sysconf(_SC_PAGESIZE);
  // Default to typical page size = 4 KB
  if(page_size <= 0)
    page_size = 4096;

  // on 64-bit systems:
  constexpr size_t page_table_entry_size = 8;

  return num_ranks * shmem_bytes * page_table_entry_size / page_size;
}
