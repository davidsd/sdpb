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

// syrk gradeschool
// IMPORTANT STEP: please make sure that a_double_array = b_double_array ahead
// of time
void mpmat::gradeschool_syrk(const int &a_start, const int &c_start,
                        const int &c_max, CBLAS_ORDER Layout,
                        CBLAS_TRANSPOSE trans, const int &n, const int &k,
                        const double alpha, const double beta) {
  int start = 2 * a_start;
  int diff = c_max - start;
  if (diff <= 2)
    {
      // base case, just multiply
      cblas_dsyrk(CblasColMajor, CblasUpper,
                  (Layout == CblasRowMajor) != (trans == CblasTrans)
                  ? CblasNoTrans
                  : CblasTrans,
                  n, k, alpha, a_double_array + k * n * a_start,
                  (Layout == CblasRowMajor) != (trans == CblasTrans) ? n : k,
                  beta, c_double_array + n * n * c_start, n);
    }
  else
    {
      // if we don't need all four, then just do grade school
    int len2 = pow(2, ceil(log2(diff)) - 1) + .1;
    int clen2 = len2 - 1; // the "fundamental length" of c one level below (post squishing)

    // C_0 = A_0 * A_0 // karatsuba
    karatsuba_syrk(a_start, c_start, c_max, Layout, trans, n, k);

    // C_1 = A_1 * A_0 // grade school
    gradeschool_gemm(a_start + len2, a_start, c_start + (3 * clen2 + 2), c_max,
                Layout, trans == CblasTrans ? CblasNoTrans : CblasTrans, trans,
                n, n, k);

    // move over C_1, deal with the overlap
    cblas_daxpy(n * n * (2 * clen2 + 1), 1.0,
                c_double_array + n * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + n * n * (c_start + clen2 + 1), 1);

    // C_1 += (A_1 * A_0)T = A_0 * A_1
    for (int i = 0; i < 2 * clen2 + 1; ++i) // in place transpose
#ifdef HAVE_MKL_H
      mkl_dimatcopy('r', 't',
#else
      cblas_dimatcopy(
          CblasRowMajor, CblasTrans,
#endif
                    n, n, 1.0,
                    c_double_array + n * n * (c_start + (3 * clen2 + 2) + i), n,
                    n);
    cblas_daxpy(n * n * (2 * clen2 + 1), 1.0,
                c_double_array + n * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + n * n * (c_start + clen2 + 1), 1);

    // clean up
    memset(c_double_array + n * n * (c_start + (3 * clen2 + 2)), 0,
           n * n * (2 * clen2 + 1) * sizeof(double));
  }
}
