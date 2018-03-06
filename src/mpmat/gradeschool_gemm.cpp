//
// Created by Rajeev Erramilli on 1/6/18.
//
#ifdef __SDPB_CUDA__
#define CUDA_API_PER_THREAD_DEFAULT_STREAM
#include "cublas_v2.h"
#include <cublasXt.h>
#include <cuda_runtime.h>
#endif

#include "mpmat.h"
#include <cassert>
#include <gmpxx.h>
#include <iostream>
#include <math.h>
//#include "tests/mpmat_tests.h"
// #include <mkl.h>
#include "../Timers.h"
#include <limits.h>
#ifdef HAVE_OMP_H
#include <omp.h>
#endif

void mpmat::gradeschool_gemm(const int &a_start, const int &b_start,
                        const int &c_start, const int &c_max,
                        CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                        CBLAS_TRANSPOSE transb, const int &m, const int &n,
                        const int &k, const double alpha, const double beta) {
  std::cout << "gradeschool gemm: " << a_start << " " << b_start << " "
            << c_start << " "
            << c_max
            << "\n";
  int start = a_start + b_start;
  int diff = c_max - start;
  // std::cerr << "gradeschool (m,n,k)=(" << m << "," << n << "," << k << ") "
  // << start << " " << c_start << " " << c_max << "\n";
  if (diff < 2) { // base case, just multiply
    // std::cerr << "attempt at gemm\n";
    // std::cerr << "\tmultiplying a[" << a_start << "] * b[" << b_start << "] =
    // c[" << c_start << "]\n";
    cblas_dgemm(Layout, transa, transb, m, n, k, alpha,
                a_double_array + k * m * a_start,
                ((Layout == CblasRowMajor) != (transa == CblasTrans)) ? k : m,
                b_double_array + k * n * b_start,
                ((Layout == CblasRowMajor) != (transb == CblasTrans)) ? n : k,
                beta, c_double_array + m * n * c_start,
                Layout == CblasRowMajor ? n : m);
  } else { // if we don't need all four, then just do grade school
    int len2 = pow(2, ceil(log2(diff)) - 1) + .1;
    int clen2 =
        len2 -
        1; // the "fundamental length" of c one level below (post squishing)
    // C_0 = A_0 * B_0 // karatsuba
    karatsuba_gemm(a_start, b_start, c_start, c_max, Layout, transa, transb, m, n,
              k);

    // C_1 = A_1 * B_0 // grade school
    gradeschool_gemm(a_start + len2, b_start, c_start + (3 * clen2 + 2), c_max,
                Layout, transa, transb, m, n, k);

    // move over C_1, deal with the overlap
    cblas_daxpy(m * n * (2 * clen2 + 1), 1.0,
                c_double_array + m * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + m * n * (c_start + clen2 + 1), 1);

    // clean up
    memset(c_double_array + m * n * (c_start + (3 * clen2 + 2)), 0,
           m * n * (2 * clen2 + 1) * sizeof(double));

    // C_2 = A_0 * B_1 // grade school // stored temporarily in C_2!
    gradeschool_gemm(a_start, b_start + len2, c_start + (3 * clen2 + 2), c_max,
                Layout, transa, transb, m, n, k);

    // C_1 += C_2
    cblas_daxpy(m * n * (2 * clen2 + 1), 1.0,
                c_double_array + m * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + m * n * (c_start + clen2 + 1), 1);

    // C_2 = 0
    memset(c_double_array + m * n * (c_start + (3 * clen2 + 2)), 0,
           m * n * (2 * clen2 + 1) * sizeof(double));
  }
}
