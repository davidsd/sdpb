//
// Created by Petr Kravchuk on 8/10/17.
//
// *****************
// *    Warning    *
// *****************
//
// We use a lot of GMP internals, and thus may be incompatible with future versions
//    Last tested version : GMP 6.1.2
//


#ifndef MPMAT_MPMAT_H
#define MPMAT_MPMAT_H

#include <gmpxx.h>
#include <mkl.h>
#include <iostream>
#include <omp.h>

#ifdef __SDPB_CUDA__
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include <cublasXt.h>
#endif

// The floating point type used to emulate integer arithmetic
typedef double mpmat_double;

// The effective number of bits in a mpmat_double mantissa
// The actual number of bits (explicitly stored) for mantissa can be lower
// We assume that mpmat_double can represent MPMAT_DOUBLE_MANT_IMPLICIT-bit integer value without loss of information
// These bit counts do not include sign bit
#define MPMAT_DOUBLE_MANT_IMPLICIT 53

// Converts a single mpf_class into an array of mpmat_double
//
// Arguments:
// (in)  source    : the mpf_class to convert
// (out) dest      : pointer to the start of the destination; must be allocated
// (in)  size      : size of the space allocated at dest, in units of mpmat_double
// (in)  mpmat_limb: the limb size to use inside mpmat_double; must be <= MPMAT_DOUBLE_MANT_IMPLICIT
// (in)  exp       : the exponent for the resulting array
//
// * dest[0] will contain the most significant limb of the number (this is reversed relative to GMP)
// * exp is the exponent, in bits, of the resulting number
// * dest is filled completely, either filling the least significant bits with 0 or discarding those of GMP,
//   depending on the relation between size and _mp_size
// * mpmat_limb should be no greater than the number of bits in mp_limb_t
// * we assume that mpmat_double 0 is represented by a zero bitstring (i.e. IEEE754)
//

class mpmat{

 private:
  //Workspace arrays for mpmat. The d_* arrays are actually arrays of pointers to the respective arrays on each GPU.
  mpmat_double *a_double_array, *b_double_array, *c_double_array, *tmp, **d_a, **d_b, **d_c;
  //this is just to keep track and avoid a segfault
  int len_a, len_b, len_c, len_t, cpu_count, gpu_count;

  #ifdef __SDPB_CUDA__
  //cuBLAS needs handles for its parameters.
  //`handle' is just there for legacy code, should
  // be removed
  cublasHandle_t handle;
  cublasHandle_t *handles;
  #endif

 public: 
  mpmat(){
    mpmat(1);
  }
 mpmat(int l){
    len_a = 0;
    len_b = 0;
    len_c = 0;
    len_t = 0;
    cpu_count = omp_get_num_threads();

    
    #ifdef __SDPB_CUDA__
    cudaGetDeviceCount(&gpu_count);
    d_a = new mpmat_double*[gpu_count];
    d_b = new mpmat_double*[gpu_count];
    d_c = new mpmat_double*[gpu_count];

    realloc_gpu(l,l,l);
    handles = new cublasHandle_t[gpu_count];

    #pragma omp parallel for
    for (int i = 0; i < gpu_count; ++i){
      cudaSetDevice(i);
      cublasCreate(handles+i);
    }
    handle = handles[0];
    #else
    realloc(l,l,l);
    gpu_count = 0;
    #endif
  }
  #ifdef __SDPB_CUDA__
  ~mpmat(){
    if (len_a != 0){
      std::cout << "deallocing a_double_array, length " << len_a << "\n";
      cudaFreeHost(a_double_array);
      for (int i = 0; i < gpu_count; ++i){
	cudaSetDevice(i);
	cudaFree(d_a[i]);
      }
    }
    if (len_b != 0){
      std::cout << "deallocing b_double_array, length " << len_b << "\n";
      cudaFreeHost(b_double_array);
      for (int i = 0; i < gpu_count; ++i){
	cudaSetDevice(i);
	cudaFree(d_b[i]);
      }
    }
    if (len_c != 0){
      std::cout << "deallocing c_double_array, length " << len_c << "\n";
      cudaFreeHost(c_double_array);
      for (int i = 0; i < gpu_count; ++i){
	cudaSetDevice(i);
	cudaFree(d_c[i]);
      }
    }
    if (len_t != 0){
      std::cout << "deallocing tmp, length " << len_t << "\n";
      delete [] tmp;
    }
    //#pragma omp parallel for
    for (int i = 0; i < gpu_count; ++i)
      cublasDestroy(handles[i]);
    delete [] handles;
    delete [] d_a;
    delete [] d_b;
    delete [] d_c;
    
  }
  void realloc_gpu(int mem_a, int mem_b, int mem_c);
  #else

~mpmat(){
    if (len_a != 0){
      std::cout << "deallocing a_double_array, length " << len_a << "\n";
      delete [] a_double_array;
    }
    if (len_b != 0){
      std::cout << "deallocing b_double_array, length " << len_b << "\n";
      delete [] b_double_array;
    }
    if (len_c != 0){
      std::cout << "deallocing c_double_array, length " << len_c << "\n";
      delete [] c_double_array;
    }
    if (len_t != 0){
      std::cout << "deallocing tmp, length " << len_t << "\n";
      delete [] tmp;
}
   
  }
  
  #endif
  
  void realloc(int mem_a, int mem_b, int mem_c);
void mpmatConvertGMPToDouble(const mpf_class source,
                             mpmat_double * dest,
                             const int size,
                             const int mpmat_limb,
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
// * mpmat_double in source can contain more bits then mpmat_limb, up to MPMAT_DOUBLE_MANT_IMPLICIT, they are extracted
//   by static_cast<mp_limb_t>; mpmat_limb essentially gives the relative exponent of different mpmat_double's
// * We assume that MPMAT_DOUBLE_MANT_IMPLICIT < mp_bits_per_limb and that mpmat_limb > 0 (strict inequalities)
// * Currently only works if dest has sufficiently many limbs for all bits in source
//     This includes + mpmat_limb * size bits for the mpmat_limbs
//                   + 1+MPMAT_DOUBLE_MANT_IMPLICIT-mpmat_limb bits for the possible overflow in the leading limb
//                   + <64 bits of padding due to exponent in GMP being proportional to 64bits
//
void mpmatConvertDoubleToGMP(mpf_class & dest,
                             const mpmat_double * source,
                             const int size,
                             const int mpmat_limb,
                             const int exp);

// Performs multilplication of GMP numbers a and b by converting to mpmat_double arrays, and stores result in dest.
// Doesn't loose precision except as described in mpmatConvertDoubleToGMP
// Test function before going to matrices
//
// Arguments:
// (out) dest : initialized mpf_class to store the result
// (in)  a    : first operand
// (in)  b    : second operand
//
void mpmatMultiplyGMPBaseCase(mpf_class & dest,
                              const mpf_class  a,
                              const mpf_class  b);

// Converts GMP vector to array of mpmat_double vectors; all operations are performed in-place over dest
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
// * This function automatically detects the maximal exponent in source and sets exp accordingly
// * This function internally uses mpmatConvertGMPToDouble
//
void mpmatConvertGMPToDoubleVector(const mpf_class * source,
                                   const int source_len,
                                   mpmat_double * dest,
                                   const int mpmat_size,
                                   const int mpmat_limb,
                                   int & exp,
                                   mpmat_double * tmp);

// Converts an array of mpmat_double vectors into a GMP vector; in-place transforms the mpmat_double array
//
// Arguments:
// (out)    dest       : the destination GMP vector; must be allocated & mpf's intialized
// (in)     dest_len   : the length of the vector
// (in)     source     : the array of source mpmat_double vectors
// (in)     mpmat_size : size of mpmat_double arrays
// (in)     mpmat_limb : the original limb size of mpmat_double
// (in)     exp        : the exponent to use for interpretation of source
// (in/out) tmp        : workplace the size of source
//
// * Note that source is transposed in-place after the call
// * This function internally uses mpmatConvertDoubleToGMP
//
void mpmatConvertDoubleToGMPVector(mpf_class * dest,
                                   const int dest_len,
                                   mpmat_double * source,
                                   const int mpmat_size,
                                   const int mpmat_limb,
                                   int exp,
                                   mpmat_double * tmp);

void mpmatConvertDoubleToGMPSymm(mpf_class * dest,
                                   const int dest_dim,
                                   mpmat_double * source,
                                   const int mpmat_size,
                                   const int mpmat_limb,
                                   int exp,
                                   mpmat_double * tmp);


void gemm_reduced(
        const CBLAS_LAYOUT Layout,
        const CBLAS_TRANSPOSE transa,
        const CBLAS_TRANSPOSE transb,
        const int m,
        const int n,
        const int k,
        //const mpf_class alpha,
        const mpf_class * a,
        //const int lda,
        const mpf_class * b,
        //const int ldb,
        //const mpf_class beta,
        mpf_class * c
        //const int ldc
);

void gemm_reduced_gpu(
        const CBLAS_LAYOUT Layout,
        const CBLAS_TRANSPOSE transa,
        const CBLAS_TRANSPOSE transb,
        const int m,
        const int n,
        const int k,
        //const mpf_class alpha,
        const mpf_class * a,
        //const int lda,
        const mpf_class * b,
        //const int ldb,
        //const mpf_class beta,
        mpf_class * c
        //const int ldc
);
void syrk_reduced_gpu(
        const CBLAS_LAYOUT Layout,
        const CBLAS_TRANSPOSE transa,
        const int m,
        const int k,
        //const mpf_class alpha,
        const mpf_class * a,
        //const int lda,
        //const mpf_class beta,
        mpf_class * c
        //const int ldc
);

void syrk_reduced(
        const CBLAS_LAYOUT Layout,
        const CBLAS_TRANSPOSE transa,
        const int m,
        const int k,
        //const mpf_class alpha,
        const mpf_class * a,
        //const int lda,
        mpf_class * c
        //const int ldc
		  );


};
 mpf_class * randomGMPVector(int size, int prec);

#endif //MPMAT_MPMAT_H
