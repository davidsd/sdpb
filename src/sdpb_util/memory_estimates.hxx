#pragma once

#include "Environment.hxx"
#include "Verbosity.hxx"
#include "assert.hxx"

#include <algorithm>

size_t
get_max_shared_memory_bytes(size_t nonshared_memory_required_per_node_bytes,
                            const Environment &env, Verbosity verbosity);

size_t bigfloat_bytes();
size_t bigfloat_bytes(mp_bitcnt_t precision);

size_t get_heap_allocated_bytes(const El::BigFloat &f);
size_t get_heap_allocated_bytes(const El::AbstractDistMatrix<El::BigFloat> &m);
size_t get_heap_allocated_bytes(const El::Matrix<El::BigFloat> &m);
template <class T> size_t get_heap_allocated_bytes(const std::vector<T> &vec)
{
  size_t res = 0;
  for(const auto &element : vec)
    res += sizeof(element) + get_heap_allocated_bytes(element);
  res += sizeof(T) * (vec.capacity() - vec.size());
  return res;
}
template <class T, std::size_t N>
size_t get_heap_allocated_bytes(const std::array<T, N> &arr)
{
  size_t res = 0;
  for(const auto &element : arr)
    res += get_heap_allocated_bytes(element);
  return res;
}
template <class T> size_t get_allocated_bytes(const T &value)
{
  return sizeof(value) + get_heap_allocated_bytes(value);
}

// Sum bytes for all ranks on a node and print from rank=0
void print_allocation_message_per_node(const Environment &env,
                                       const std::string &name, size_t bytes);

// Extra memory required to compute L^{-1} X, via
// El::Trsm(LEFT, LOWER, NORMAL, NON_UNIT, 1, L, X);
// L ~ height x height, lower triangular
// X ~ height x width
// See https://gitlab.com/bootstrapcollaboration/elemental/-/blob/3ee9f404e632d9d0a4751e9b6717da51e59b5697/src/blas_like/level3/Trsm.cpp#L119-L122
//
// NB: if grid_size (=MPI group size) is a prime number (e.g. 11, 13, or 17),
// then the grid has 1D shape (e.g. 11x1, 13x1, or 17x1).
// In this case the terms proportional to grid_height will lead to
// very high memory consumption compared to more balanced grids (e.g. 4x3 or 4x4).
// To avoid this issue, set e.g. --procGranularity 2 or --procGranularity 4.
inline size_t get_trsm_bytes(const int height, const int width,
                             const int grid_height, const int grid_width)
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

// Extra memory required to compute X L, via
// El::Trmm(RIGHT, LOWER, NORMAL, NON_UNIT, 1, L, X)
// X ~ height x width
// L ~ width x width, lower triangular
// See https://gitlab.com/bootstrapcollaboration/elemental/-/blob/3ee9f404e632d9d0a4751e9b6717da51e59b5697/src/blas_like/level3/Trmm.cpp#L78
//
// NB: if grid_size (=MPI group size) is a prime number (e.g. 11, 13, or 17),
// then the grid has 1D shape (e.g. 11x1, 13x1, or 17x1).
// In this case the terms proportional to grid_height will lead to
// very high memory consumption compared to more balanced grids (e.g. 4x3 or 4x4).
// To avoid this issue, set e.g. --procGranularity 2 or --procGranularity 4.
inline size_t get_trmm_bytes(const int height, const int width,
                             const int grid_height, const int grid_width)
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
