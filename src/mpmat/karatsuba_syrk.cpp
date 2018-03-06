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

void mpmat::karatsuba_syrk(const int &a_start, const int &c_start, const int &c_max,
                      CBLAS_ORDER Layout, CBLAS_TRANSPOSE trans, const int &n,
                      const int &k, const double alpha, const double beta) {

  int start = 2 * a_start;
  int diff = c_max - start;
  std::cout << "karatsuba syrk: " << a_start << " " << c_start << " "
            << c_max << " " << diff
            << "\n";
  
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
      // if we need all four quadrants, do full karatsuba
    int len2 = pow(2, ceil(log2(diff)) - 2) + .1; // the "fundamental length" of a and b
    int clen2 = len2 - 1; // the "fundamental length" of c one level below (post squishing)

    // C_0 = A_0 * A_0 // full karatsuba
    karatsuba_syrk(a_start, c_start, 2 * len2 + start, Layout, trans, n, k);

    std::cout << "karatsuba syrk A_O+=A_1: " << a_start << " " << c_start << " "
              << c_max << "\n";
    // A_0 += A_1
    cblas_daxpy(k * n * len2, 1.0, a_double_array + k * n * (a_start + len2), 1,
                a_double_array + k * n * a_start, 1);

    // C_1 = A_0 * A_0 // full karatsuba
    karatsuba_syrk(a_start, c_start + (2 * clen2 + 1), 2 * len2 + start, Layout,
              trans, n, k);

    std::cout << "karatsuba syrk A_O-=A_1: " << a_start << " " << c_start << " "
              << c_max << "\n";
    // A_0 -= A_1 // clean up
    cblas_daxpy(k * n * len2, -1.0, a_double_array + k * n * (a_start + len2),
                1, a_double_array + k * n * a_start, 1);

    // C_1 -= C_0
    cblas_daxpy(n * n * (2 * clen2 + 1), -1.0, c_double_array + n * n * c_start,
                1, c_double_array + n * n * (c_start + (2 * clen2 + 1)), 1);

    // C_02 += C_10
    cblas_daxpy(n * n * clen2, 1.0,
                c_double_array + n * n * (c_start + 2 * clen2 + 1), 1,
                c_double_array + n * n * (c_start + clen2 + 1), 1);

    // squish the rest of C_1 against C_0
    cblas_dcopy(n * n * clen2,
                c_double_array + n * n * (c_start + 3 * clen2 + 1), 1,
                c_double_array + n * n * (c_start + 2 * clen2 + 1), 1);
    cblas_dcopy(n * n, c_double_array + n * n * (c_start + 4 * clen2 + 1), 1,
                c_double_array + n * n * (c_start + 3 * clen2 + 1), 1);

    // clear out the old
    memset(c_double_array + n * n * (c_start + 3 * clen2 + 2), 0,
           n * n * clen2 * sizeof(double));

    // C_2 = A_1 * A_1 // conditional karatsuba, depending on the cutoff
    if (diff <=
        3 * len2) // if we need half or less of C_2, then just use grade school
      gradeschool_syrk(a_start + len2, c_start + (3 * clen2 + 2), c_max, Layout,
                  trans, n, k);
    else // otherwise, we need all quadrants and therefore karatsuba
      karatsuba_syrk(a_start + len2, c_start + (3 * clen2 + 2), c_max, Layout, trans,
                n, k);

    std::cout << "karatsuba syrk C_1-=C_2: " << a_start << " " << c_start << " "
              << c_max << "\n";
    // C_1 -= C_2
    cblas_daxpy(n * n * (2 * clen2 + 1), -1.0,
                c_double_array + n * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + n * n * (c_start + clen2 + 1), 1);

    // C_12 += C_20
    cblas_daxpy(n * n * clen2, 1.0,
                c_double_array + n * n * (c_start + 3 * clen2 + 2), 1,
                c_double_array + n * n * (c_start + 2 * clen2 + 2), 1);

    // squish the rest
    cblas_dcopy(n * n * clen2,
                c_double_array + n * n * (c_start + 4 * clen2 + 2), 1,
                c_double_array + n * n * (c_start + 3 * clen2 + 2), 1);
    cblas_dcopy(n * n, c_double_array + n * n * (c_start + 5 * clen2 + 2), 1,
                c_double_array + n * n * (c_start + 4 * clen2 + 2), 1);

    memset(c_double_array + n * n * (c_start + 4 * clen2 + 3), 0,
           n * n * clen2 * sizeof(double));
  }
}
