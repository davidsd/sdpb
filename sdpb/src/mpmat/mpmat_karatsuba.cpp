//
// Created by Rajeev Erramilli on 1/6/18.
//
#ifdef __SDPB_CUDA__
#define CUDA_API_PER_THREAD_DEFAULT_STREAM
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include <cublasXt.h>
#endif

#include "mpmat.h"
#include <gmpxx.h>
#include <math.h>
#include <cassert>
#include <iostream>
//#include "tests/mpmat_tests.h"
#include <mkl.h>
#include <limits.h>
#include "../Timers.h"
#include <omp.h>

template <typename T>
inline T min(T a, T b) { return a < b ? a : b ; }

template <typename T>
inline T max (T a,T b) { return a>b ? a : b; }

//  max(mem_t,pow(3,ceil(log2(mem_c))));

void mpmat::karatsuba(int a_start, int b_start, int c_start, int c_max,
                      CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE transa,
                      CBLAS_TRANSPOSE transb, int m, int n, int k,
                      const double alpha = 1.0, const double beta = 1.0){
  int start = a_start + b_start;
  int len = pow(2,ceil(log2(c_max - start))-1); //review
  if(len == 0){ // base case, just multiply
    cblas_dgemm(
            Layout,
            transa,
            transb,
            m,
            n,
            k,
            alpha,
            a_double_array+k*m*a_start,
            Layout == CblasRowMajor ? k : m,
            b_double_array+k*n*b_start,
            Layout == CblasRowMajor ? n : k,
            beta,
            c_double_array+m*n*c_start,
            Layout == CblasRowMajor ? n : m
    );
  }
  else if (c_max > 2*len + start){ //if we need all four quadrants, do full karatsuba
    // C_0 = A_0 * B_0 // full karatsuba
    karatsuba(a_start, b_start, c_start, 2*len + start, Layout, transa, transb, m, n, k);
    // C_2 = A_1 * B_1 // conditional karatsuba
    karatsuba(a_start+len, b_start+len, c_start+2*len, c_max, Layout, transa, transb, m, n, k);
    // A_0 += A_1
    cblas_daxpy(k*m, 1.0, a_double_array+k*m*(a_start+len), 1, a_double_array+k*m*a_start, 1);
    // B_0 += B_1
    cblas_daxpy(k*n, 1.0, b_double_array+k*m*(b_start+len), 1, b_double_array+k*m*b_start, 1);
    // C_1 = A_0 * B_0 // full karatsuba
    karatsuba(a_start, b_start, c_start+len, 2*len + start, Layout, transa, transb, m, n, k);
    // C_1 -= C_0
    cblas_daxpy(m*n, -1.0, c_double_array+m*n*c_start, 1, c_double_array+m*n*(c_start+len), 1);
    // C_1 -= C_2
    cblas_daxpy(m*n, -1.0, c_double_array+m*n*(c_start+2*len), 1, c_double_array+m*n*(c_start+len), 1);
    // A_0 -= A_1 // clean up
    cblas_daxpy(k*m, -1.0, a_double_array+k*m*a_start, 1, a_double_array+k*m*(a_start+len), 1);
    // B_0 -= B_1 // clean up
    cblas_daxpy(k*n, -1.0, b_double_array+k*n*b_start, 1, b_double_array+k*n*(b_start+len), 1);
  }
  else { //if we don't need all four, then just do grade school
    // C_0 = A_0 * B_0 // conditional karatsuba
    karatsuba(a_start, b_start, c_start, c_max, Layout, transa, transb, m, n, k);
    // C_1 = A_1 * B_0 // conditional karatsuba
    karatsuba(a_start+len, b_start, c_start+len, c_max, Layout, transa, transb, m, n, k);
    // C_2 = A_0 * B_1 // conditional karatsuba // stored in C_2 temporarily!
    karatsuba(a_start, b_start+len, c_start+2*len, c_max, Layout, transa, transb, m, n, k);
    // C_1 += C_2
    cblas_daxpy(m*n, 1.0, c_double_array+m*n*(c_start+2*len), 1, c_double_array+m*n*(c_start+len), 1);
    // C_2 = 0
    memset(c_double_array+m*n*(c_start+2*len),0,m*n*len*sizeof(double));
  }
}

void mpmat::karatsuba(int a_start, int c_start, int c_max,
                      CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE trans, int n, int k,
                      const double alpha = 1.0, const double beta = 1.0){
  int start = 2*a_start;
  int len = pow(2,ceil(log2(c_max - start))-1); //review. Make sure we're getting all cases accurately.
  if(len == 0){ // base case, just multiply
    cblas_dsyrk(CblasColMajor,
                CblasUpper,(Layout == CblasRowMajor) != (transa == CblasTrans) ? CblasNoTrans : CblasTrans,
                n,k,alpha,
                a_double_array+k*n*a_start,
                (Layout == CblasRowMajor) != (transa == CblasTrans) ? m : k,
                beta,
                c_double_array+n*n*c_start,
                n);
  }
  else if (c_max > 2*len + start){ //if we need all four quadrants, do full symmetric karatsuba
    // C_0 = A_0 * A_0 // full karatsuba
    karatsuba(a_start, c_start, 2*len + start, Layout, trans, n, k);
    // C_2 = A_1 * A_1 // conditional karatsuba
    karatsuba(a_start+len, c_start+2*len, c_max, Layout, trans, n, k);
    // A_0 += A_1
    cblas_daxpy(k*n, 1.0, a_double_array+k*n*(a_start+len), 1, a_double_array+k*n*a_start, 1);
    // C_1 = A_0 * A_0 // full karatsuba
    karatsuba(a_start, c_start+len, 2*len + start, Layout, trans, n, k);
    // C_1 -= C_0
    cblas_daxpy(m*n, -1.0, c_double_array+n*n*c_start, 1, c_double_array+n*n*(c_start+len), 1);
    // C_1 -= C_2
    cblas_daxpy(m*n, -1.0, c_double_array+n*n*(c_start+2*len), 1, c_double_array+n*n*(c_start+len), 1);
    // A_0 -= A_1 // clean up
    cblas_daxpy(k*n, -1.0, a_double_array+k*n*a_start, 1, a_double_array+k*n*(a_start+len), 1);
  }
  else { //if we don't need all four, then just do symmetric grade school
    // C_0 = A_0 * B_0 // conditional karatsuba
    karatsuba(a_start, c_max, Layout, transa, transb, m, n, k);
    // C_1 = A_1 * B_0 // conditional karatsuba
    karatsuba(a_start+len, b_start, c_start+len, c_max, Layout, transa, transb, m, n, k);
    // C_1 += C_2 // NEED TO TRANSPOSE
    cblas_daxpy(m*n, 1.0, c_double_array+m*n*(c_start+2*len), 1, c_double_array+m*n*(c_start+len), 1);
    // C_2 = 0
    memset(c_double_array+m*n*(c_start+2*len),0,m*n*len*sizeof(double));
  }
}

