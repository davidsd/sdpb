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
#include <omp.h>

void mpmat::karatsuba_gemm(const int &a_start, const int &b_start,
                      const int &c_start, const int &c_max, CBLAS_ORDER Layout,
                      CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb,
                      const int &m, const int &n, const int &k,
                      const double alpha, const double beta) {
  std::cout << "karatsuba gemm: " << a_start << " " << b_start << " "
            << c_start << " "
            << c_max
            << "\n";
  int start = a_start + b_start;
  int diff = c_max - start;
  if (diff <= 2) { // base case, just multiply
    cblas_dgemm(Layout, transa, transb, m, n, k, alpha,
                a_double_array + k * m * a_start,
                ((Layout == CblasRowMajor) != (transa == CblasTrans)) ? k : m,
                b_double_array + k * n * b_start,
                ((Layout == CblasRowMajor) != (transb == CblasTrans)) ? n : k,
                beta, c_double_array + m * n * c_start,
                Layout == CblasRowMajor ? n : m);
  } else { // if we need all four quadrants, do full karatsuba
    int len2 = pow(2, ceil(log2(diff)) - 2) + .1; // the "fundamental length" of a and b
    int clen2 = len2 - 1; // the "fundamental length" of c one level below (post squishing)

    // C_0 = A_0 * B_0 // full karatsuba
    karatsuba_gemm(a_start, b_start, c_start, 2 * len2 + start, Layout, transa,
              transb, m, n, k);

    // A_0 += A_1
    cblas_daxpy(k * m * len2, 1.0, a_double_array + k * m * (a_start + len2), 1,
                a_double_array + k * m * a_start, 1);

    // B_0 += B_1
    cblas_daxpy(k * n * len2, 1.0, b_double_array + k * n * (b_start + len2), 1,
                b_double_array + k * n * b_start, 1);

    // C_1 = A_0 * B_0 // full karatsuba
    karatsuba_gemm(a_start, b_start, c_start + (2 * clen2 + 1), 2 * len2 + start,
              Layout, transa, transb, m, n, k);

    // A_0 -= A_1 // clean up
    cblas_daxpy(k * m * len2, -1.0, a_double_array + k * m * (a_start + len2),
                1, a_double_array + k * m * a_start, 1);

    // B_0 -= B_1 // clean up
    cblas_daxpy(k * n * len2, -1.0, b_double_array + k * n * (b_start + len2),
                1, b_double_array + k * n * b_start, 1);

    // C_1 -= C_0
    cblas_daxpy(m * n * (2 * clen2 + 1), -1.0, c_double_array + m * n * c_start,
                1, c_double_array + m * n * (c_start + (2 * clen2 + 1)), 1);

    // C_02 += C_10
    cblas_daxpy(m * n * clen2, 1.0,
                c_double_array + m * n * (c_start + 2 * clen2 + 1), 1,
                c_double_array + m * n * (c_start + clen2 + 1), 1);

    // squish the rest of C_1 against C_0
    cblas_dcopy(m * n * clen2,
                c_double_array + m * n * (c_start + 3 * clen2 + 1), 1,
                c_double_array + m * n * (c_start + 2 * clen2 + 1), 1);
    cblas_dcopy(m * n, c_double_array + m * n * (c_start + 4 * clen2 + 1), 1,
                c_double_array + m * n * (c_start + 3 * clen2 + 1), 1);

    // clear out the old
    memset(c_double_array + m * n * (c_start + 3 * clen2 + 2), 0,
           m * n * clen2 * sizeof(double));

    // C_2 = A_1 * B_1 // conditional karatsuba, depending on the cutoff
    if (diff <=
        3 * len2) // if we need half or less of C_2, then just use grade school
      gradeschool_gemm(a_start + len2, b_start + len2, c_start + (3 * clen2 + 2),
                  c_max, Layout, transa, transb, m, n, k);
    else // otherwise, we need all quadrants and therefore karatsuba
      karatsuba_gemm(a_start + len2, b_start + len2, c_start + (3 * clen2 + 2),
                c_max, Layout, transa, transb, m, n, k);

    // C_1 -= C_2
    cblas_daxpy(m * n * (2 * clen2 + 1), -1.0,
                c_double_array + m * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + m * n * (c_start + clen2 + 1), 1);

    // C_12 += C_20
    cblas_daxpy(m * n * clen2, 1.0,
                c_double_array + m * n * (c_start + 3 * clen2 + 2), 1,
                c_double_array + m * n * (c_start + 2 * clen2 + 2), 1);

    // squish the rest
    cblas_dcopy(m * n * clen2,
                c_double_array + m * n * (c_start + 4 * clen2 + 2), 1,
                c_double_array + m * n * (c_start + 3 * clen2 + 2), 1);
    cblas_dcopy(m * n, c_double_array + m * n * (c_start + 5 * clen2 + 2), 1,
                c_double_array + m * n * (c_start + 4 * clen2 + 2), 1);

    memset(c_double_array + m * n * (c_start + 4 * clen2 + 3), 0,
           m * n * clen2 * sizeof(double));
  }
}
