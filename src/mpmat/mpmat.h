//
// Created by Petr Kravchuk on 8/10/17.
//
// *****************
// *    Warning    *
// *****************
//
// We use a lot of GMP internals, and thus may be incompatible with future
// versions
//    Last tested version : GMP 6.1.2
//

#ifndef MPMAT_MPMAT_H
#define MPMAT_MPMAT_H

#include <gmpxx.h>

#ifdef HAVE_MKL_H
#include <mkl.h>
#else
#include <cblas.h>
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

#ifdef __SDPB_CUDA__
#include "cublas_v2.h"
#include <cuda_runtime.h>
#endif

// The floating point type used to emulate integer arithmetic
typedef double mpmat_double;

// The effective number of bits in a mpmat_double mantissa
// The actual number of bits (explicitly stored) for mantissa can be lower
// We assume that mpmat_double can represent MPMAT_DOUBLE_MANT_IMPLICIT-bit
// integer value without loss of information These bit counts do not include
// sign bit
#define MPMAT_DOUBLE_MANT_IMPLICIT 53

// Converts a single mpf_class into an array of mpmat_double
//
// Arguments:
// (in)  source    : the mpf_class to convert
// (out) dest      : pointer to the start of the destination; must be allocated
// (in)  size      : size of the space allocated at dest, in units of
// mpmat_double (in)  mpmat_limb: the limb size to use inside mpmat_double; must
// be <= MPMAT_DOUBLE_MANT_IMPLICIT (in)  exp       : the exponent for the
// resulting array
//
// * dest[0] will contain the most significant limb of the number (this is
// reversed relative to GMP)
// * exp is the exponent, in bits, of the resulting number
// * dest is filled completely, either filling the least significant bits with 0
// or discarding those of GMP,
//   depending on the relation between size and _mp_size
// * mpmat_limb should be no greater than the number of bits in mp_limb_t
// * we assume that mpmat_double 0 is represented by a zero bitstring (i.e.
// IEEE754)
//

class mpmat {

private:
  // Workspace arrays for mpmat. The d_* arrays are actually arrays of pointers
  // to the respective arrays on each GPU.

  // this is just to keep track and avoid a segfault
  size_t len_a, len_b, len_c, len_t, *gpu_len_a, *gpu_len_b, *gpu_len_c;
  int cpu_count, gpu_count;

#ifdef __SDPB_CUDA__
  // cuBLAS needs handles for its parameters.
  //`handle' is just there for legacy code, should
  // be removed
  size_t gpu_mem;
  cublasHandle_t handle;
  cublasHandle_t *handles;
#endif

#ifdef __SDPB_CUDA__
  void karatsuba_gpu(const int &a_start, const int &b_start, const int &c_start,
                     const int &c_max, CBLAS_ORDER Layout,
                     CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb,
                     const int &m, const int &n, const int &k, int device = 0,
                     const double alpha = 1.0, const double beta = 1.0);
  // Implements a symmetric version of the above
  void karatsuba_gpu(const int &a_start, const int &c_start, const int &c_max,
                     CBLAS_ORDER Layout, CBLAS_TRANSPOSE trans, const int &n,
                     const int &k, int device = 0, const double alpha = 1.0,
                     const double beta = 1.0);
  void karatsuba_cpu(const int &a_start, const int &b_start, const int &c_start,
                     const int &c_max, CBLAS_ORDER Layout,
                     CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb,
                     const int &m, const int &n, const int &k,
                     const double alpha = 1.0, const double beta = 1.0);
  // Implements a symmetric version of the above
  void karatsuba_cpu(const int &a_start, const int &c_start, const int &c_max,
                     CBLAS_ORDER Layout, CBLAS_TRANSPOSE trans, const int &n,
                     const int &k, const double alpha = 1.0,
                     const double beta = 1.0);
#endif
  void karatsuba_gemm(const int &a_start, const int &b_start, const int &c_start,
                 const int &c_max, CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                 CBLAS_TRANSPOSE transb, const int &m, const int &n,
                 const int &k, const double alpha = 1.0,
                 const double beta = 1.0);
  // Implements a symmetric version of the above
  void karatsuba_syrk(const int &a_start, const int &c_start, const int &c_max,
                 CBLAS_ORDER Layout, CBLAS_TRANSPOSE trans, const int &n,
                 const int &k, const double alpha = 1.0,
                 const double beta = 1.0);

  void karatsuba_bc(const int &a_start, const int &b_start, const int &c_start,
                    const int &c_max, CBLAS_ORDER Layout,
                    CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb,
                    const int &m, const int &n, const int &k,
                    const double alpha = 1.0, const double beta = 1.0);

#ifdef __SDPB_CUDA__
  // Implements a truncated recursive multiplication that maximizes Karatsuba
  void gradeschool_gpu(const int &a_start, const int &b_start,
                       const int &c_start, const int &c_max, CBLAS_ORDER Layout,
                       CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb,
                       const int &m, const int &n, const int &k, int device = 0,
                       const double alpha = 1.0, const double beta = 1.0);
  // Implements a symmetric version of the above
  void gradeschool_gpu(const int &a_start, const int &c_start, const int &c_max,
                       CBLAS_ORDER Layout, CBLAS_TRANSPOSE trans, const int &n,
                       const int &k, int device = 0, const double alpha = 1.0,
                       const double beta = 1.0);
  void gradeschool_cpu(const int &a_start, const int &b_start,
                       const int &c_start, const int &c_max, CBLAS_ORDER Layout,
                       CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb,
                       const int &m, const int &n, const int &k,
                       const double alpha = 1.0, const double beta = 1.0);
  // Implements a symmetric version of the above
  void gradeschool_cpu(const int &a_start, const int &c_start, const int &c_max,
                       CBLAS_ORDER Layout, CBLAS_TRANSPOSE trans, const int &n,
                       const int &k, const double alpha = 1.0,
                       const double beta = 1.0);
#endif
  // Implements a truncated recursive multiplication that maximizes Karatsuba
  void gradeschool_gemm(const int &a_start, const int &b_start, const int &c_start,
                   const int &c_max, CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                   CBLAS_TRANSPOSE transb, const int &m, const int &n,
                   const int &k, const double alpha = 1.0,
                   const double beta = 1.0);
  // Implements a symmetric version of the above
  void gradeschool_syrk(const int &a_start, const int &c_start, const int &c_max,
                   CBLAS_ORDER Layout, CBLAS_TRANSPOSE trans, const int &n,
                   const int &k, const double alpha = 1.0,
                   const double beta = 1.0);

  void gradeschool_bc(const int &a_start, const int &b_start,
                      const int &c_start, const int &c_max, CBLAS_ORDER Layout,
                      CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb,
                      const int &m, const int &n, const int &k,
                      const double alpha = 1.0, const double beta = 1.0);

  void clear_gpu(int device);

public:
  mpmat_double *a_double_array, *b_double_array, *c_double_array, *tmp, **d_a,
      **d_b, **d_c;
  mpmat(int l = 1) {
    len_a = 0;
    len_b = 0;
    len_c = 0;
    len_t = 0;
#ifdef HAVE_OMP_H
    cpu_count = omp_get_num_threads();
#else
    cpu_count = 1;
#endif
    
#ifdef __SDPB_CUDA__
    cudaGetDeviceCount(&gpu_count);
    d_a = new mpmat_double *[gpu_count];
    d_b = new mpmat_double *[gpu_count];
    d_c = new mpmat_double *[gpu_count];
    gpu_len_a = new size_t[gpu_count];
    gpu_len_b = new size_t[gpu_count];
    gpu_len_c = new size_t[gpu_count];

    realloc_gpu(l, l, l);
    handles = new cublasHandle_t[gpu_count];

    size_t total;

    for (int i = 0; i < gpu_count; ++i) {
      cudaSetDevice(i);
      cublasCreate(handles + i);
      size_t tmp_free, tmp_total;
      cudaMemGetInfo(&tmp_free, &tmp_total);
      if (i == 0)
        total = tmp_total;
      else
        assert(tmp_total == total);
      gpu_mem = tmp_free; // DANGEROUS: one GPU might randomly have
                          // significantly less GPU memory free.
    }
    gpu_mem *= .75; // safety buffer
    std::cerr << "there is " << gpu_mem / (1 << 20)
              << " MiB of gpu memory available.\n";
    handle = handles[0];
#else
    realloc(l, l, l);
    gpu_count = 0;
#endif
  }
#ifdef __SDPB_CUDA__
  ~mpmat() {
    if (len_a != 0) {
      cudaFreeHost(a_double_array);
      for (int i = 0; i < gpu_count; ++i) {
        cudaSetDevice(i);
        cudaFree(d_a[i]);
      }
    }
    if (len_b != 0) {
      cudaFreeHost(b_double_array);
      for (int i = 0; i < gpu_count; ++i) {
        cudaSetDevice(i);
        cudaFree(d_b[i]);
      }
    }
    if (len_c != 0) {
      cudaFreeHost(c_double_array);
      for (int i = 0; i < gpu_count; ++i) {
        cudaSetDevice(i);
        cudaFree(d_c[i]);
      }
    }
    if (len_t != 0) {
      delete[] tmp;
    }
    //#pragma omp parallel for
    for (int i = 0; i < gpu_count; ++i)
      cublasDestroy(handles[i]);
    delete[] handles;
    delete[] d_a;
    delete[] d_b;
    delete[] d_c;
    delete[] gpu_len_a;
    delete[] gpu_len_b;
    delete[] gpu_len_c;
  }
  void realloc_gpu(size_t mem_a, size_t mem_b, size_t mem_c);
  void realloc_gpu_only(size_t mem_a, size_t mem_b, size_t mem_c,
                        int device = 0);
#else

  ~mpmat() {
    if (len_a != 0) {
      delete[] a_double_array;
    }
    if (len_b != 0) {
      delete[] b_double_array;
    }
    if (len_c != 0) {
      delete[] c_double_array;
    }
    if (len_t != 0) {
      delete[] tmp;
    }
  }

#endif

  void realloc(size_t mem_a, size_t mem_b, size_t mem_c);
  void mpmatConvertGMPToDouble(const mpf_class source, mpmat_double *dest,
                               const int size, const int mpmat_limb,
                               const int exp);

  // Converts an array of mpmat_double into mpf_class
  //
  // Arguments:
  // (out) dest       : the destination mpf_variable
  // (in)  source     : pointer the start of the mpmat_double array
  // (in)  size       : length of source
  // (in)  mpmat_limb : the original limb size of the source
  // (in)  exp        : the exponent of the array, in bits
  //
  // * mpmat_double in source can contain more bits then mpmat_limb, up to
  // MPMAT_DOUBLE_MANT_IMPLICIT, they are extracted
  //   by static_cast<mp_limb_t>; mpmat_limb essentially gives the relative
  //   exponent of different mpmat_double's
  // * We assume that MPMAT_DOUBLE_MANT_IMPLICIT < mp_bits_per_limb and that
  // mpmat_limb > 0 (strict inequalities)
  // * Currently only works if dest has sufficiently many limbs for all bits in
  // source
  //     This includes + mpmat_limb * size bits for the mpmat_limbs
  //                   + 1+MPMAT_DOUBLE_MANT_IMPLICIT-mpmat_limb bits for the
  //                   possible overflow in the leading limb
  //                   + <64 bits of padding due to exponent in GMP being
  //                   proportional to 64bits
  //
  void mpmatConvertDoubleToGMP(mpf_class &dest, const mpmat_double *source,
                               const int size, const int mpmat_limb,
                               const int exp);

  // Performs multilplication of GMP numbers a and b by converting to
  // mpmat_double arrays, and stores result in dest. Doesn't loose precision
  // except as described in mpmatConvertDoubleToGMP Test function before going
  // to matrices
  //
  // Arguments:
  // (out) dest : initialized mpf_class to store the result
  // (in)  a    : first operand
  // (in)  b    : second operand
  //
  void mpmatMultiplyGMPBaseCase(mpf_class &dest, const mpf_class a,
                                const mpf_class b);

  // Converts GMP vector to array of mpmat_double vectors; all operations are
  // performed in-place over dest
  //
  // Arguments:
  // (in)  source     : the pointer to the first element in the GMP vector
  // (in)  source_len : the length of the GMP vector
  // (out) dest       : the pointer to the allocated array of mpmat_doubles
  // (in)  size       : number of mpmat_doubles to use per GMP entry
  // (in)  mpmat_limb : the bit limb size to use in mpmat_double array
  // (out) exp        : the exponent used in the output
  // (in/out) tmp     : workplace the size of dest
  //
  // * The same exponent exp is used for all entries of the vector
  // * This function automatically detects the maximal exponent in source and
  // sets exp accordingly
  // * This function internally uses mpmatConvertGMPToDouble
  //
  void mpmatConvertGMPToDoubleVector(const mpf_class *source,
                                     const int source_len, mpmat_double *dest,
                                     const int mpmat_size, const int mpmat_limb,
                                     int &exp, mpmat_double *tmp);

  // Converts an array of mpmat_double vectors into a GMP vector; in-place
  // transforms the mpmat_double array
  //
  // Arguments:
  // (out)    dest       : the destination GMP vector; must be allocated & mpf's
  // intialized (in)     dest_len   : the length of the vector (in)     source
  // : the array of source mpmat_double vectors (in)     mpmat_size : size of
  // mpmat_double arrays (in)     mpmat_limb : the original limb size of
  // mpmat_double (in)     exp        : the exponent to use for interpretation
  // of source (in/out) tmp        : workplace the size of source
  //
  // * Note that source is transposed in-place after the call
  // * This function internally uses mpmatConvertDoubleToGMP
  //
  void mpmatConvertDoubleToGMPVector(mpf_class *dest, const int dest_len,
                                     mpmat_double *source, const int mpmat_size,
                                     const int mpmat_limb, int exp,
                                     mpmat_double *tmp);

  void mpmatConvertDoubleToGMPSymm(mpf_class *dest, const int dest_dim,
                                   mpmat_double *source, const int mpmat_size,
                                   const int mpmat_limb, int exp,
                                   mpmat_double *tmp);

  void gemm_reduced(const CBLAS_ORDER Layout, const CBLAS_TRANSPOSE transa,
                    const CBLAS_TRANSPOSE transb, const int m, const int n,
                    const int k,
                    // const mpf_class alpha,
                    const mpf_class *a,
                    // const int lda,
                    const mpf_class *b,
                    // const int ldb,
                    // const mpf_class beta,
                    mpf_class *c
                    // const int ldc
  );

  void gemm_reduced_gpu(const CBLAS_ORDER Layout, const CBLAS_TRANSPOSE transa,
                        const CBLAS_TRANSPOSE transb, const int m, const int n,
                        const int k,
                        // const mpf_class alpha,
                        const mpf_class *a,
                        // const int lda,
                        const mpf_class *b,
                        // const int ldb,
                        // const mpf_class beta,
                        mpf_class *c
                        // const int ldc
  );
  void syrk_reduced_gpu(const CBLAS_ORDER Layout, const CBLAS_TRANSPOSE transa,
                        const int m, const int k,
                        // const mpf_class alpha,
                        const mpf_class *a,
                        // const int lda,
                        // const mpf_class beta,
                        mpf_class *c
                        // const int ldc
  );

  void syrk_reduced(const CBLAS_ORDER Layout, const CBLAS_TRANSPOSE transa,
                    const size_t m, const size_t k,
                    // const mpf_class alpha,
                    const mpf_class *a,
                    // const int lda,
                    mpf_class *c
                    // const int ldc
  );
// Implements a recursive Karatsuba multiplication on an array of matrices with
// a cutoff c_max
#ifdef __SDPB_CUDA__
  void karatsuba_generic(const int &c_max, CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                 CBLAS_TRANSPOSE transb, const int &m, const int &n,
                 const int &k, bool gpu = true);
  // symmetric case
  void karatsuba_symmetric(const int &c_max, CBLAS_ORDER Layout, CBLAS_TRANSPOSE trans,
                 const int &n, const int &k, bool gpu = true);
#else
  void karatsuba_generic(const int &c_max, CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                 CBLAS_TRANSPOSE transb, const int &m, const int &n,
                 const int &k);
  // symmetric case
  void karatsuba_symmetric(const int &c_max, CBLAS_ORDER Layout, CBLAS_TRANSPOSE trans,
                 const int &n, const int &k);
#endif

  void treecondense(double *c, int size, int l);

  void mpmat_conversion_test(int i, int f, int d);

  bool karatsuba_test(int m, int n, int k, int l);
  bool symm_karatsuba_test(int n, int k, int l);
#ifdef __SDPB_CUDA__
  bool karatsuba_test_gpu(int m, int n, int k, int l);
  bool symm_karatsuba_test_gpu(int n, int k, int l);
#endif
  bool base_karatsuba_test();
};

bool compareSymmMatrices(double *x, double *y, int n, int l, bool lower = true,
                         bool rowmaj = true);
mpf_class *randomGMPVector(int size, int prec);
void print_mpf_bits(const mpf_class &a);
bool compare_mpf_bits(const mpf_class &a, const mpf_class &b);

void print_mpmat_double_array(const mpmat_double *array, int len);
#endif // MPMAT_MPMAT_H
