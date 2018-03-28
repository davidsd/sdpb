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

template <typename T> T min(T a, T b) { return a < b ? a : b; }

template <typename T> T max(T a, T b) { return a > b ? a : b; }

//  max(mem_t,pow(3,ceil(log2(mem_c))));
#ifdef __SDPB_CUDA__
void mpmat::clear_gpu(int device) {
  cudaMemset(d_a[device], 0, gpu_len_a[device] * sizeof(mpmat_double));
  cudaMemset(d_b[device], 0, gpu_len_b[device] * sizeof(mpmat_double));
  cudaMemset(d_c[device], 0, gpu_len_c[device] * sizeof(mpmat_double));
}

void mpmat::karatsuba(const int &c_max, CBLAS_ORDER Layout,
                      CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb,
                      const int &m, const int &n, const int &k, bool gpu) {
  if (gpu) {
    gradeschool_cpu(0, 0, 0, c_max, Layout, transa, transb, m, n, k);
  }
  // gradeschool_gpu(0, 0, 0, c_max, Layout, transa, transb, m, n, k);
  else
    gradeschool(0, 0, 0, c_max, Layout, transa, transb, m, n, k);
}
#else
void mpmat::karatsuba_generic(const int &c_max, CBLAS_ORDER Layout,
                      CBLAS_TRANSPOSE transa, CBLAS_TRANSPOSE transb,
                      const int &m, const int &n, const int &k) {
  gradeschool_gemm(0, 0, 0, c_max, Layout, transa, transb, m, n, k);
}
#endif

#ifdef __SDPB_CUDA__
void mpmat::karatsuba(const int &c_max, CBLAS_ORDER Layout,
                      CBLAS_TRANSPOSE trans, const int &n, const int &k,
                      bool gpu) {
  if (gpu) {
    double *b_tmp = b_double_array;
    b_double_array = a_double_array;
    gradeschool_cpu(0, 0, c_max, Layout, trans, n, k);
    // std::cerr << "really done with recursings\n";
    b_double_array = b_tmp;
  }
  // double * b_tmp[gpu_count];
  //   for (int i = 0; i < gpu_count; ++i){
  //     b_tmp[i] = d_b[i];
  //     d_b[i] = d_a[i];
  //   }
  //   gradeschool_gpu(0, 0, c_max, Layout, trans, n, k);
  //   //std::cerr << "really done with recursings\n";
  //   for (int i = 0; i < gpu_count; ++i)
  //     d_b[i] = b_tmp[i];
  else {
    double *b_tmp = b_double_array;
    b_double_array = a_double_array;
    gradeschool(0, 0, c_max, Layout, trans, n, k);
    // std::cerr << "really done with recursings\n";
    b_double_array = b_tmp;
  }
}
#else
void mpmat::karatsuba_symmetric(const int &c_max, CBLAS_ORDER Layout,
                      CBLAS_TRANSPOSE trans, const int &n, const int &k) {
  double *b_tmp = b_double_array;
  b_double_array = a_double_array;
  gradeschool_syrk(0, 0, c_max, Layout, trans, n, k);
  // std::cerr << "really done with recursings\n";
  b_double_array = b_tmp;
}
#endif

#ifdef __SDPB_CUDA__
// gemm karatsuba
void mpmat::karatsuba_cpu(const int &a_start, const int &b_start,
                          const int &c_start, const int &c_max,
                          CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                          CBLAS_TRANSPOSE transb, const int &m, const int &n,
                          const int &k, const double alpha, const double beta) {
  int start = a_start + b_start;
  int diff = c_max - start;
  int len2 =
      pow(2, ceil(log2(diff)) - 2) + .1; // the "fundamental length" of a and b
  int clen2 =
      len2 -
      1; // the "fundamental length" of c one level below (post squishing)
  size_t mem = (len2 * m * k + len2 * n * k + (6 * len2 + 2) * m * n) *
               sizeof(mpmat_double);
  if (diff <= 2) { // base case, just multiply. Currently set to do so on CPU:
                   // investigate switching to GPU in some case?
    cblas_dgemm(Layout, transa, transb, m, n, k, alpha,
                a_double_array + k * m * a_start,
                ((Layout == CblasRowMajor) != (transa == CblasTrans)) ? k : m,
                b_double_array + k * n * b_start,
                ((Layout == CblasRowMajor) != (transb == CblasTrans)) ? n : k,
                beta, c_double_array + m * n * c_start,
                Layout == CblasRowMajor ? n : m);
  } else if (mem < gpu_mem) { // if the problem can fit on the GPU in the next
                              // step, do so.

    for (int i = 0; i < 2; ++i) {
      realloc_gpu_only(len2 * m * k, len2 * n * k, (6 * len2 + 2) * m * n, i);
      clear_gpu(i);
      cudaMemcpy(d_a[i], a_double_array + k * m * (a_start + i * len2),
                 len2 * m * k * sizeof(mpmat_double), cudaMemcpyHostToDevice);
      cudaMemcpy(d_b[i], b_double_array + k * n * (b_start + i * len2),
                 len2 * n * k * sizeof(mpmat_double), cudaMemcpyHostToDevice);
    }
    // A_0 += A_1
    cblas_daxpy(k * m * len2, 1.0, a_double_array + k * m * (a_start + len2), 1,
                a_double_array + k * m * a_start, 1);

    // B_0 += B_1
    cblas_daxpy(k * n * len2, 1.0, b_double_array + k * n * (b_start + len2), 1,
                b_double_array + k * n * b_start, 1);
    {
      { // C_0 = A_0 * B_0 // full karatsuba
        cudaSetDevice(0);
        karatsuba_gpu(0, 0, 0, 2 * len2, Layout, transa, transb, m, n, k, 0);
      }
      {
        cudaSetDevice(1);
        if (diff <= 3 * len2) { // if we need half or less of C_2, then just use
                                // grade school
          std::cerr << "gs ";
          gradeschool_gpu(0, 0, 0, diff - 2 * len2, Layout, transa, transb, m,
                          n, k, 1); // Maybe?
        } else { // otherwise, we need all quadrants and therefore karatsuba
          std::cerr << "ka ";
          karatsuba_gpu(0, 0, 0, diff - 2 * len2, Layout, transa, transb, m, n,
                        k, 1); // Maybe?
        }
      }
      {
        realloc_gpu_only(len2 * m * k, len2 * n * k, (6 * len2 + 2) * m * n, 2);
        cudaMemcpy(d_a[2], a_double_array + k * m * (a_start),
                   len2 * m * k * sizeof(mpmat_double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_b[2], b_double_array + k * n * (b_start),
                   len2 * n * k * sizeof(mpmat_double), cudaMemcpyHostToDevice);

        // C_1 = A_0 * B_0 // full karatsuba
        karatsuba_gpu(0, 0, 0, 2 * len2, Layout, transa, transb, m, n, k, 2);
      }
    }
    // A_0 -= A_1 // clean up
    cblas_daxpy(k * m * len2, -1.0, a_double_array + k * m * (a_start + len2),
                1, a_double_array + k * m * a_start, 1);

    // B_0 -= B_1 // clean up
    cblas_daxpy(k * n * len2, -1.0, b_double_array + k * n * (b_start + len2),
                1, b_double_array + k * n * b_start, 1);

    cudaSetDevice(0);
    cudaMemcpy(c_double_array + m * n * c_start, d_c[0],
               (2 * clen2 + 1) * m * n * sizeof(mpmat_double),
               cudaMemcpyDeviceToHost);
    cudaSetDevice(1);
    cudaMemcpy(c_double_array + m * n * (c_start + 2 * clen2 + 1), d_c[1],
               (2 * clen2 + 1) * m * n * sizeof(mpmat_double),
               cudaMemcpyDeviceToHost);

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
    cudaSetDevice(2);
    cudaMemcpy(c_double_array + m * n * (c_start + 3 * clen2 + 2), d_c[2],
               (2 * clen2 + 1) * m * n * sizeof(mpmat_double),
               cudaMemcpyDeviceToHost);

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
  } else { // if there's not enough room, keep recursing

    // C_0 = A_0 * B_0 // full karatsuba
    karatsuba_cpu(a_start, b_start, c_start, 2 * len2 + start, Layout, transa,
                  transb, m, n, k);

    // A_0 += A_1
    cblas_daxpy(k * m * len2, 1.0, a_double_array + k * m * (a_start + len2), 1,
                a_double_array + k * m * a_start, 1);

    // B_0 += B_1
    cblas_daxpy(k * n * len2, 1.0, b_double_array + k * n * (b_start + len2), 1,
                b_double_array + k * n * b_start, 1);

    // C_1 = A_0 * B_0 // full karatsuba
    karatsuba_cpu(a_start, b_start, c_start + (2 * clen2 + 1), 2 * len2 + start,
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
      gradeschool_cpu(a_start + len2, b_start + len2, c_start + (3 * clen2 + 2),
                      c_max, Layout, transa, transb, m, n, k);
    else // otherwise, we need all quadrants and therefore karatsuba
      karatsuba_cpu(a_start + len2, b_start + len2, c_start + (3 * clen2 + 2),
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
// gemm karatsuba
void mpmat::karatsuba_gpu(const int &a_start, const int &b_start,
                          const int &c_start, const int &c_max,
                          CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                          CBLAS_TRANSPOSE transb, const int &m, const int &n,
                          const int &k, int device, const double alpha,
                          const double beta) {
  cudaSetDevice(device);
  // std::cerr << "karatsuba " << a_start << "," << b_start << "," << c_start <<
  // "," << c_max << "," << m << "," << n << "," << k << "," << device << "\n";
  // std::cerr << (((transa == CblasTrans) != (Layout == CblasColMajor)) ?
  // "CUBLAS_OP_N " : "CUBLAS_OP_T " )<< (((transb == CblasTrans) != (Layout ==
  // CblasColMajor)) ? "CUBLAS_OP_N\n" : "CUBLAS_OP_T\n"); std::cerr <<
  // (((Layout == CblasColMajor) != (transa == CblasTrans)) ? k : n) << " " <<
  // (((Layout == CblasColMajor) != (transb == CblasTrans)) ? m : k )<< "\n";
  int start = a_start + b_start;
  int diff = c_max - start;
  if (diff <= 2) { // base case, just multiply
    // std::cerr << "basecase'd out\n";
    cublasDgemm(
        handles[device],
        ((transb == CblasTrans) != (Layout == CblasColMajor)) ? CUBLAS_OP_T
                                                              : CUBLAS_OP_N,
        ((transa == CblasTrans) != (Layout == CblasColMajor)) ? CUBLAS_OP_T
                                                              : CUBLAS_OP_N,
        m, n, k, &alpha, d_b[device] + k * n * b_start,
        ((Layout == CblasColMajor) != (transb == CblasTrans)) ? k : n,
        d_a[device] + k * m * a_start,
        ((Layout == CblasColMajor) != (transa == CblasTrans)) ? m : k, &beta,
        d_c[device] + m * n * c_start, Layout == CblasColMajor ? n : m);
  } else { // if we need all four quadrants, do full karatsuba
    int len2 = pow(2, ceil(log2(diff)) - 2) +
               .1; // the "fundamental length" of a and b
    int clen2 =
        len2 -
        1; // the "fundamental length" of c one level below (post squishing)
    double minus = -1.0;
    // std::cerr << "len2 = " << len2 << "\n";
    // C_0 = A_0 * B_0 // full karatsuba
    karatsuba_gpu(a_start, b_start, c_start, 2 * len2 + start, Layout, transa,
                  transb, m, n, k, device);

    // A_0 += A_1
    cublasDaxpy(handles[device], k * m * len2, &alpha,
                d_a[device] + k * m * (a_start + len2), 1,
                d_a[device] + k * m * a_start, 1);

    // B_0 += B_1
    cublasDaxpy(handles[device], k * n * len2, &alpha,
                d_b[device] + k * n * (b_start + len2), 1,
                d_b[device] + k * n * b_start, 1);

    // C_1 = A_0 * B_0 // full karatsuba
    karatsuba_gpu(a_start, b_start, c_start + (2 * clen2 + 1), 2 * len2 + start,
                  Layout, transa, transb, m, n, k, device);

    // A_0 -= A_1 // clean up
    cublasDaxpy(handles[device], k * m * len2, &minus,
                d_a[device] + k * m * (a_start + len2), 1,
                d_a[device] + k * m * a_start, 1);

    // B_0 -= B_1 // clean up
    cublasDaxpy(handles[device], k * n * len2, &minus,
                d_b[device] + k * n * (b_start + len2), 1,
                d_b[device] + k * n * b_start, 1);

    // C_1 -= C_0
    cublasDaxpy(handles[device], m * n * (2 * clen2 + 1), &minus,
                d_c[device] + m * n * c_start, 1,
                d_c[device] + m * n * (c_start + (2 * clen2 + 1)), 1);

    // C_02 += C_10
    cublasDaxpy(handles[device], m * n * clen2, &alpha,
                d_c[device] + m * n * (c_start + 2 * clen2 + 1), 1,
                d_c[device] + m * n * (c_start + clen2 + 1), 1);

    // squish the rest of C_1 against C_0
    cublasDcopy(handles[device], m * n * clen2,
                d_c[device] + m * n * (c_start + 3 * clen2 + 1), 1,
                d_c[device] + m * n * (c_start + 2 * clen2 + 1), 1);
    cublasDcopy(handles[device], m * n,
                d_c[device] + m * n * (c_start + 4 * clen2 + 1), 1,
                d_c[device] + m * n * (c_start + 3 * clen2 + 1), 1);

    // clear out the old
    cudaMemset(d_c[device] + m * n * (c_start + 3 * clen2 + 2), 0,
               m * n * clen2 * sizeof(double));

    // C_2 = A_1 * B_1 // conditional karatsuba, depending on the cutoff
    if (diff <=
        3 * len2) // if we need half or less of C_2, then just use grade school
      gradeschool_gpu(a_start + len2, b_start + len2, c_start + (3 * clen2 + 2),
                      c_max, Layout, transa, transb, m, n, k, device);
    else // otherwise, we need all quadrants and therefore karatsuba
      karatsuba_gpu(a_start + len2, b_start + len2, c_start + (3 * clen2 + 2),
                    c_max, Layout, transa, transb, m, n, k, device);

    // C_1 -= C_2
    cublasDaxpy(handles[device], m * n * (2 * clen2 + 1), &minus,
                d_c[device] + m * n * (c_start + (3 * clen2 + 2)), 1,
                d_c[device] + m * n * (c_start + clen2 + 1), 1);

    // C_12 += C_20
    cublasDaxpy(handles[device], m * n * clen2, &alpha,
                d_c[device] + m * n * (c_start + 3 * clen2 + 2), 1,
                d_c[device] + m * n * (c_start + 2 * clen2 + 2), 1);

    // squish the rest
    cublasDcopy(handles[device], m * n * clen2,
                d_c[device] + m * n * (c_start + 4 * clen2 + 2), 1,
                d_c[device] + m * n * (c_start + 3 * clen2 + 2), 1);
    cublasDcopy(handles[device], m * n,
                d_c[device] + m * n * (c_start + 5 * clen2 + 2), 1,
                d_c[device] + m * n * (c_start + 4 * clen2 + 2), 1);

    cudaMemset(d_c[device] + m * n * (c_start + 4 * clen2 + 3), 0,
               m * n * clen2 * sizeof(double));
  }
}
#endif

#ifdef __SDPB_CUDA__
void mpmat::karatsuba_cpu(const int &a_start, const int &c_start,
                          const int &c_max, CBLAS_ORDER Layout,
                          CBLAS_TRANSPOSE trans, const int &n, const int &k,
                          const double alpha, const double beta) {

  int start = 2 * a_start;
  int diff = c_max - start;
  int len2 =
      pow(2, ceil(log2(diff)) - 2) + .1; // the "fundamental length" of a and b
  int clen2 =
      len2 -
      1; // the "fundamental length" of c one level below (post squishing)
  size_t mem = (len2 * n * k + (6 * len2 + 2) * n * n) * sizeof(mpmat_double);
  if (diff <= 2) { // base case, just multiply
    cblas_dsyrk(CblasColMajor, CblasUpper,
                (Layout == CblasRowMajor) != (trans == CblasTrans)
                    ? CblasNoTrans
                    : CblasTrans,
                n, k, alpha, a_double_array + k * n * a_start,
                (Layout == CblasRowMajor) != (trans == CblasTrans) ? n : k,
                beta, c_double_array + n * n * c_start, n);
  } else if (mem < gpu_mem) {

    for (int i = 0; i < 3; ++i) {
      realloc_gpu_only(len2 * n * k, 0, (6 * len2 + 2) * n * n, i);
      clear_gpu(i);
      cudaMemcpy(d_a[i], a_double_array + k * n * (a_start + i * len2),
                 len2 * n * k * sizeof(mpmat_double), cudaMemcpyHostToDevice);
    }
    // A_0 += A_1
    cblas_daxpy(k * n * len2, 1.0, a_double_array + k * n * (a_start + len2), 1,
                a_double_array + k * n * a_start, 1);

    {
      { // C_0 = A_0 * B_0 // full karatsuba
        cudaSetDevice(0);
        double *d_tmp_b = d_b[0];
        d_b[0] = d_a[0];
        karatsuba_gpu(0, 0, 2 * len2, Layout, trans, n, k, 0);
        d_b[0] = d_tmp_b;
      }
      {
        cudaSetDevice(1);
        double *d_tmp_b = d_b[1];
        d_b[1] = d_a[1];
        if (diff <= 3 * len2) { // if we need half or less of C_2, then just use
                                // grade school
          std::cerr << "gs ";
          gradeschool_gpu(0, 0, diff - 2 * len2, Layout, trans, n, k,
                          1); // Maybe?
        } else { // otherwise, we need all quadrants and therefore karatsuba
          std::cerr << "ka ";
          karatsuba_gpu(0, 0, diff - 2 * len2, Layout, trans, n, k, 1); // Maybe?
        }
        d_b[1] = d_tmp_b;
      }
      {
        realloc_gpu_only(len2 * n * k, 0, (6 * len2 + 2) * n * n, 2);
        cudaMemcpy(d_a[2], a_double_array + k * n * (a_start),
                   len2 * n * k * sizeof(mpmat_double), cudaMemcpyHostToDevice);

        double *d_tmp_b = d_b[2];
        d_b[2] = d_a[2];
        // C_1 = A_0 * B_0 // full karatsuba
        karatsuba_gpu(0, 0, 2 * len2, Layout, trans, n, k, 2);
        d_b[2] = d_tmp_b;
      }
    }
    // A_0 -= A_1 // clean up
    cblas_daxpy(k * n * len2, -1.0, a_double_array + k * n * (a_start + len2),
                1, a_double_array + k * n * a_start, 1);

    // C_0 = A_0 * A_0 // full karatsuba
    cudaSetDevice(0);
    cudaMemcpy(c_double_array + n * n * c_start, d_c[0],
               (2 * clen2 + 1) * n * n * sizeof(mpmat_double),
               cudaMemcpyDeviceToHost);
    // C_1 = A_0 * A_0 // full karatsuba
    cudaSetDevice(1);
    cudaMemcpy(c_double_array + n * n * (c_start + 2 * clen2 + 1), d_c[1],
               (2 * clen2 + 1) * n * n * sizeof(mpmat_double),
               cudaMemcpyDeviceToHost);

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
    cudaSetDevice(2);
    cudaMemcpy(c_double_array + n * n * (c_start + 3 * clen2 + 2), d_c[2],
               (2 * clen2 + 1) * n * n * sizeof(mpmat_double),
               cudaMemcpyDeviceToHost);

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
  } else { // if we need all four quadrants, do full karatsuba
    // C_0 = A_0 * A_0 // full karatsuba
    karatsuba_cpu(a_start, c_start, 2 * len2 + start, Layout, trans, n, k);

    // A_0 += A_1
    cblas_daxpy(k * n * len2, 1.0, a_double_array + k * n * (a_start + len2), 1,
                a_double_array + k * n * a_start, 1);

    // C_1 = A_0 * A_0 // full karatsuba
    karatsuba_cpu(a_start, c_start + (2 * clen2 + 1), 2 * len2 + start, Layout,
                  trans, n, k);

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
      gradeschool_cpu(a_start + len2, c_start + (3 * clen2 + 2), c_max, Layout,
                      trans, n, k);
    else // otherwise, we need all quadrants and therefore karatsuba
      karatsuba_cpu(a_start + len2, c_start + (3 * clen2 + 2), c_max, Layout,
                    trans, n, k);

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
// syrk karatsuba
void mpmat::karatsuba_gpu(const int &a_start, const int &c_start,
                          const int &c_max, CBLAS_ORDER Layout,
                          CBLAS_TRANSPOSE trans, const int &n, const int &k,
                          int device, const double alpha, const double beta) {
  // std::cerr << "karatsuba on the " << device << " gpu.\n";
  // std::cerr << "skaratsuba " << a_start << "," << c_start << "," << c_max <<
  // "," << n << "," << k << "\n";
  int start = 2 * a_start;
  int diff = c_max - start;
  cudaSetDevice(device);
  if (diff <= 2) { // base case, just multiply
    cublasDsyrk(handles[device], CUBLAS_FILL_MODE_UPPER,
                (Layout == CblasRowMajor) != (trans == CblasTrans)
                    ? CUBLAS_OP_N
                    : CUBLAS_OP_T,
                n, k, &alpha, d_a[device] + k * n * a_start,
                (Layout == CblasRowMajor) != (trans == CblasTrans) ? n : k,
                &beta, d_c[device] + n * n * c_start, n);
    // if (stat == CUBLAS_STATUS_SUCCESS) std::cerr << "multiplication
    // succeeded\n"; else if (stat == CUBLAS_STATUS_NOT_INITIALIZED) std::cerr
    // << "cublas wasn't set up\n"; else std::cerr << stat << "\n";

  } else { // if we need all four quadrants, do full karatsuba
    int len2 = pow(2, ceil(log2(diff)) - 2) +
               .1; // the "fundamental length" of a and b
    int clen2 =
        len2 -
        1; // the "fundamental length" of c one level below (post squishing)
    double minus = -1.0;

    // C_0 = A_0 * A_0 // full karatsuba
    karatsuba_gpu(a_start, c_start, 2 * len2 + start, Layout, trans, n, k,
                  device);

    // A_0 += A_1
    cublasDaxpy(handles[device], k * n * len2, &alpha,
                d_a[device] + k * n * (a_start + len2), 1,
                d_a[device] + k * n * a_start, 1);

    // C_1 = A_0 * A_0 // full karatsuba
    karatsuba_gpu(a_start, c_start + (2 * clen2 + 1), 2 * len2 + start, Layout,
                  trans, n, k, device);

    // A_0 -= A_1 // clean up
    cublasDaxpy(handles[device], k * n * len2, &minus,
                d_a[device] + k * n * (a_start + len2), 1,
                d_a[device] + k * n * a_start, 1);

    // C_1 -= C_0
    cublasDaxpy(handles[device], n * n * (2 * clen2 + 1), &minus,
                d_c[device] + n * n * c_start, 1,
                d_c[device] + n * n * (c_start + (2 * clen2 + 1)), 1);

    // C_02 += C_10
    cublasDaxpy(handles[device], n * n * clen2, &alpha,
                d_c[device] + n * n * (c_start + 2 * clen2 + 1), 1,
                d_c[device] + n * n * (c_start + clen2 + 1), 1);

    // squish the rest of C_1 against C_0
    cublasDcopy(handles[device], n * n * clen2,
                d_c[device] + n * n * (c_start + 3 * clen2 + 1), 1,
                d_c[device] + n * n * (c_start + 2 * clen2 + 1), 1);
    cublasDcopy(handles[device], n * n,
                d_c[device] + n * n * (c_start + 4 * clen2 + 1), 1,
                d_c[device] + n * n * (c_start + 3 * clen2 + 1), 1);

    // clear out the old
    cudaMemset(d_c[device] + n * n * (c_start + 3 * clen2 + 2), 0,
               n * n * clen2 * sizeof(double));

    // C_2 = A_1 * A_1 // conditional karatsuba, depending on the cutoff
    if (diff <=
        3 * len2) // if we need half or less of C_2, then just use grade school
      gradeschool_gpu(a_start + len2, c_start + (3 * clen2 + 2), c_max, Layout,
                      trans, n, k, device);
    else // otherwise, we need all quadrants and therefore karatsuba
      karatsuba_gpu(a_start + len2, c_start + (3 * clen2 + 2), c_max, Layout,
                    trans, n, k, device);

    // C_1 -= C_2
    cublasDaxpy(handles[device], n * n * (2 * clen2 + 1), &minus,
                d_c[device] + n * n * (c_start + (3 * clen2 + 2)), 1,
                d_c[device] + n * n * (c_start + clen2 + 1), 1);

    // C_12 += C_20
    cublasDaxpy(handles[device], n * n * clen2, &alpha,
                d_c[device] + n * n * (c_start + 3 * clen2 + 2), 1,
                d_c[device] + n * n * (c_start + 2 * clen2 + 2), 1);

    // squish the rest
    cublasDcopy(handles[device], n * n * clen2,
                d_c[device] + n * n * (c_start + 4 * clen2 + 2), 1,
                d_c[device] + n * n * (c_start + 3 * clen2 + 2), 1);
    cublasDcopy(handles[device], n * n,
                d_c[device] + n * n * (c_start + 5 * clen2 + 2), 1,
                d_c[device] + n * n * (c_start + 4 * clen2 + 2), 1);

    cudaMemset(d_c[device] + n * n * (c_start + 4 * clen2 + 3), 0,
               n * n * clen2 * sizeof(double));
  }
}
#endif
#ifdef __SDPB_CUDA__
// gemm gradeschool
void mpmat::gradeschool_cpu(const int &a_start, const int &b_start,
                            const int &c_start, const int &c_max,
                            CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                            CBLAS_TRANSPOSE transb, const int &m, const int &n,
                            const int &k, const double alpha,
                            const double beta) {
  int start = a_start + b_start;
  int diff = c_max - start;
  int len2 = pow(2, ceil(log2(diff)) - 1) + .1;
  int clen2 =
      len2 -
      1; // the "fundamental length" of c one level below (post squishing)
  size_t mem = (len2 * m * k + len2 * n * k + (6 * len2 + 2) * m * n) *
               sizeof(mpmat_double);
  // std::cerr << "gradeschool (m,n,k)=(" << m << "," << n << "," << k << ") "
  // << start << " " << c_start << " " << c_max << "\n"; std::cerr << "len2 = "
  // << len2 << "\n";
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
  } else if (mem < gpu_mem) {
    for (int i = 0; i < 3; ++i) {
      realloc_gpu_only(len2 * m * k, len2 * n * k, (6 * len2 + 2) * m * n, i);
      clear_gpu(i);
      cudaSetDevice(i);
      cudaMemcpy(d_a[i], a_double_array + k * m * (a_start + (i % 2) * len2),
                 len2 * m * k * sizeof(mpmat_double), cudaMemcpyHostToDevice);
      cudaMemcpy(d_b[i], b_double_array + k * n * (b_start + (i / 2) * len2),
                 len2 * n * k * sizeof(mpmat_double), cudaMemcpyHostToDevice);
    }
// for (int i = 0; i < 3; ++i){
//   realloc_gpu_only(len2*m*k,len2*n*k,(6*len2+2)*m*n,0);
//   cudaSetDevice(0);
//   cudaMemcpy(d_a[0],a_double_array+k*m*(a_start+(i % 2)*len2),
//   len2*m*k*sizeof(mpmat_double),cudaMemcpyHostToDevice);
//   cudaMemcpy(d_b[0],b_double_array+k*n*(b_start+(i / 2)*len2),
//   len2*n*k*sizeof(mpmat_double),cudaMemcpyHostToDevice);
// }
    {
      { // C_0 = A_0 * B_0 // full karatsuba
        cudaSetDevice(0);
        karatsuba_gpu(0, 0, 0, 2 * len2, Layout, transa, transb, m, n, k, 0);
        cudaDeviceSynchronize();
        cudaMemcpy(c_double_array + m * n * c_start, d_c[0],
                   (2 * clen2 + 1) * m * n * sizeof(mpmat_double),
                   cudaMemcpyDeviceToHost);
      }
      {
        cudaSetDevice(1);
        gradeschool_gpu(0, 0, 0, diff - len2, Layout, transa, transb, m, n, k,
                        1); // Maybe?
        // cudaMemcpy(c_double_array+m*n*(c_start+3*clen2+2), d_c[1],
        // (2*clen2+1)*m*n*sizeof(mpmat_double),cudaMemcpyDeviceToHost);
      }
      {
        cudaSetDevice(2);
        gradeschool_gpu(0, 0, 0, diff - len2, Layout, transa, transb, m, n, k,
                        2); // Maybe?
      }
    }
    cudaDeviceSynchronize();
    // cudaSetDevice(0);
    // cudaMemcpy(c_double_array+m*n*c_start, d_c[0],
    // (2*clen2+1)*m*n*sizeof(mpmat_double),cudaMemcpyDeviceToHost);
    memset(c_double_array + m * n * (c_max), 0,
           m * n * (2 * clen2 + 1) * sizeof(double));
    cudaSetDevice(1);
    cudaMemcpy(c_double_array + m * n * (c_start + 3 * clen2 + 2), d_c[1],
               (2 * clen2 + 1) * m * n * sizeof(mpmat_double),
               cudaMemcpyDeviceToHost);

    // move over C_1, deal with the overlap
    cblas_daxpy(m * n * (2 * clen2 + 1), 1.0,
                c_double_array + m * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + m * n * (c_start + clen2 + 1), 1);

    // clean up
    memset(c_double_array + m * n * (c_start + (3 * clen2 + 2)), 0,
           m * n * (2 * clen2 + 1) * sizeof(double));

    cudaSetDevice(2);
    cudaMemcpy(c_double_array + m * n * (c_start + 3 * clen2 + 2), d_c[2],
               (2 * clen2 + 1) * m * n * sizeof(mpmat_double),
               cudaMemcpyDeviceToHost);

    // C_1 += C_2
    cblas_daxpy(m * n * (2 * clen2 + 1), 1.0,
                c_double_array + m * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + m * n * (c_start + clen2 + 1), 1);

    // C_2 = 0
    memset(c_double_array + m * n * (c_start + (3 * clen2 + 2)), 0,
           m * n * (2 * clen2 + 1) * sizeof(double));
  } else { // if we don't need all four, then just do grade school

    // C_0 = A_0 * B_0 // karatsuba
    karatsuba_cpu(a_start, b_start, c_start, c_max, Layout, transa, transb, m,
                  n, k);

    // C_1 = A_1 * B_0 // grade school
    gradeschool_cpu(a_start + len2, b_start, c_start + (3 * clen2 + 2), c_max,
                    Layout, transa, transb, m, n, k);

    // move over C_1, deal with the overlap
    cblas_daxpy(m * n * (2 * clen2 + 1), 1.0,
                c_double_array + m * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + m * n * (c_start + clen2 + 1), 1);

    // clean up
    memset(c_double_array + m * n * (c_start + (3 * clen2 + 2)), 0,
           m * n * (2 * clen2 + 1) * sizeof(double));

    // C_2 = A_0 * B_1 // grade school // stored temporarily in C_2!
    gradeschool_cpu(a_start, b_start + len2, c_start + (3 * clen2 + 2), c_max,
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
// gemm gradeschool
void mpmat::gradeschool_gpu(const int &a_start, const int &b_start,
                            const int &c_start, const int &c_max,
                            CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                            CBLAS_TRANSPOSE transb, const int &m, const int &n,
                            const int &k, int device, const double alpha,
                            const double beta) {
  // std::cerr << "gradeschool " << a_start << "," << b_start << "," << c_start
  // << "," << c_max << "," << m << "," << n << "," << k << "\n"; std::cerr <<
  // (((transa == CblasTrans) != (Layout == CblasColMajor)) ? "CUBLAS_OP_N " :
  // "CUBLAS_OP_T " )<< (((transb == CblasTrans) != (Layout == CblasColMajor)) ?
  // "CUBLAS_OP_N\n" : "CUBLAS_OP_T\n"); std::cerr << (((Layout ==
  // CblasColMajor) != (transa == CblasTrans)) ? k : n) << " " << (((Layout ==
  // CblasColMajor) != (transb == CblasTrans)) ? m : k )<< "\n";
  cudaSetDevice(device);
  int start = a_start + b_start;
  int diff = c_max - start;
  if (diff < 2) { // base case, just multiply
    cublasDgemm(
        handles[device],
        ((transb == CblasTrans) != (Layout == CblasColMajor)) ? CUBLAS_OP_T
                                                              : CUBLAS_OP_N,
        ((transa == CblasTrans) != (Layout == CblasColMajor)) ? CUBLAS_OP_T
                                                              : CUBLAS_OP_N,
        m, n, k, &alpha, d_b[device] + k * n * b_start,
        ((Layout == CblasColMajor) != (transb == CblasTrans)) ? k : n,
        d_a[device] + k * m * a_start,
        ((Layout == CblasColMajor) != (transa == CblasTrans)) ? m : k, &beta,
        d_c[device] + m * n * c_start, Layout == CblasColMajor ? n : m);
  } else { // if we don't need all four, then just do grade school
    // std::cerr << "testing\n";
    int len2 = pow(2, ceil(log2(diff)) - 1) + .1;
    int clen2 =
        len2 -
        1; // the "fundamental length" of c one level below (post squishing)
    // C_0 = A_0 * B_0 // karatsuba
    karatsuba_gpu(a_start, b_start, c_start, c_max, Layout, transa, transb, m,
                  n, k, device);

    // C_1 = A_1 * B_0 // grade school
    gradeschool_gpu(a_start + len2, b_start, c_start + (3 * clen2 + 2), c_max,
                    Layout, transa, transb, m, n, k, device);

    // move over C_1, deal with the overlap
    // std::cerr << "adding c[" << c_start+(3*clen2+2) << ":" <<
    // c_start+(3*clen2+2) + (2*clen2+1) << "] to c[" << (c_start+clen2+1) << ":"
    // << (c_start+clen2+1) + (2*clen2+1) << "]\n";
    cublasDaxpy(handles[device], m * n * (2 * clen2 + 1), &alpha,
                d_c[device] + m * n * (c_start + (3 * clen2 + 2)), 1,
                d_c[device] + m * n * (c_start + clen2 + 1), 1);

    // clean up
    // std::cerr << "clearing c[" << (c_start+(3*clen2+2)) << ":" <<
    // (c_start+(3*clen2+2)) + (2*clen2+1) << "]\n";
    cudaMemset(d_c[device] + m * n * (c_start + (3 * clen2 + 2)), 0,
               m * n * (2 * clen2 + 1) * sizeof(double));

    // C_2 = A_0 * B_1 // grade school // stored temporarily in C_2!
    gradeschool_gpu(a_start, b_start + len2, c_start + (3 * clen2 + 2), c_max,
                    Layout, transa, transb, m, n, k, device);

    // C_1 += C_2
    // std::cerr << "adding c[" << c_start+(3*clen2+2) << ":" <<
    // c_start+(3*clen2+2) +(2*clen2+1) << "] to c[" << (c_start+clen2+1) << ":"
    // <<  (c_start+clen2+1)+(2*clen2+1) << "]\n";
    cublasDaxpy(handles[device], m * n * (2 * clen2 + 1), &alpha,
                d_c[device] + m * n * (c_start + (3 * clen2 + 2)), 1,
                d_c[device] + m * n * (c_start + clen2 + 1), 1);

    // C_2 = 0
    // std::cerr << "clearing c[" << (c_start+(3*clen2+2)) << ":" <<
    // (c_start+(3*clen2+2)) + (2*clen2+1) << "]\n";
    cudaMemset(d_c[device] + m * n * (c_start + (3 * clen2 + 2)), 0,
               m * n * (2 * clen2 + 1) * sizeof(double));
  }
}
#endif

#ifdef __SDPB_CUDA__
// syrk gradeschool
// IMPORTANT STEP: please make sure that a_double_array = b_double_array ahead
// of time
void mpmat::gradeschool_cpu(const int &a_start, const int &c_start,
                            const int &c_max, CBLAS_ORDER Layout,
                            CBLAS_TRANSPOSE trans, const int &n, const int &k,
                            const double alpha, const double beta) {
  int start = 2 * a_start;
  int diff = c_max - start;
  int len2 = pow(2, ceil(log2(diff)) - 1) + .1;
  int clen2 =
      len2 -
      1; // the "fundamental length" of c one level below (post squishing)
  size_t mem =
      (2 * len2 * n * k + (6 * len2 + 2) * n * n) * sizeof(mpmat_double);
  // std::cerr << "sgradeschool (n,k)=(" << n << "," << k << ") " << start << "
  // " << c_start << " " << c_max << "\n";
  if (diff < 2) { // base case, just multiply
    // std::cerr << "basecase'd out\n";
    cblas_dsyrk(CblasColMajor, CblasUpper,
                (Layout == CblasRowMajor) != (trans == CblasTrans)
                    ? CblasNoTrans
                    : CblasTrans,
                n, k, alpha, a_double_array + k * n * a_start,
                (Layout == CblasRowMajor) != (trans == CblasTrans) ? n : k,
                beta, c_double_array + n * n * c_start, n);
  } else if (mem < gpu_mem) {
    // std::cerr << "PUTTING STUFF ON THE GPUS\n";
    cudaSetDevice(3);
    realloc_gpu_only(len2 * n * k, 0, (6 * len2 + 2) * n * n, 3);
    clear_gpu(3);
    cudaMemcpy(d_a[3], a_double_array + k * n * (a_start),
               len2 * n * k * sizeof(mpmat_double), cudaMemcpyHostToDevice);
    {
      { // C_0 = A_0 * B_0 // full karatsuba
        cudaSetDevice(3);
        double *d_tmp_b = d_b[3];
        d_b[3] = d_a[3];
        karatsuba_gpu(0, 0, 2 * len2, Layout, trans, n, k, 3);
        d_b[3] = d_tmp_b;
        cudaMemcpy(c_double_array + n * n * c_start, d_c[3],
                   (2 * clen2 + 1) * n * n * sizeof(mpmat_double),
                   cudaMemcpyDeviceToHost);
        // memset(c_double_array+n*n*(c_max),0,n*n*(2*clen2+1)*sizeof(double));
      }
      { // do the rest in the next level down
        gradeschool_cpu(
            a_start + len2, a_start, c_start + (3 * clen2 + 2), c_max, Layout,
            trans == CblasTrans ? CblasNoTrans : CblasTrans, trans, n, n, k);
      }
    }
    // C_0 = A_0 * A_0 // karatsuba

    memset(c_double_array + n * n * (c_max), 0,
           n * n * (2 * clen2 + 1) * sizeof(double));
    // happened already
    // gradeschool happened already

    // move over C_1, deal with the overlap
    // std::cerr << "\tadding c[" << c_start+(3*clen2+2) << ":" <<
    // c_start+(3*clen2+2)+2*clen2+1 << "] -> c[" <<(c_start+clen2+1) << ":" <<
    // (c_start+3*clen2+2) <<"]\n";
    cblas_daxpy(n * n * (2 * clen2 + 1), 1.0,
                c_double_array + n * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + n * n * (c_start + clen2 + 1), 1);

    // C_1 += (A_1 * A_0)T = A_0 * A_1
    // std::cerr << "TRANSPOSING\n";
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
    // std::cerr << "TRANSPOSED\n";
    // std::cerr << "\tadding c[" << c_start+(3*clen2+2) << ":" <<
    // c_start+(3*clen2+2)+2*clen2+1 << "]T -> c[" <<(c_start+clen2+1) << ":" <<
    // (c_start+3*clen2+2) <<"]\n";
    cblas_daxpy(n * n * (2 * clen2 + 1), 1.0,
                c_double_array + n * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + n * n * (c_start + clen2 + 1), 1);

    // clean up
    memset(c_double_array + n * n * (c_start + (3 * clen2 + 2)), 0,
           n * n * (2 * clen2 + 1) * sizeof(double));
  } else { // if we don't need all four, then just do grade school

    // C_0 = A_0 * A_0 // karatsuba
    // std::cerr << "\tmultiplying a[" << a_start << ":" << a_start+len2 << "]
    // -> c[" << c_start << ":" << c_start+ 2*clen2 + 1 << "]\n";
    karatsuba_cpu(a_start, c_start, c_max, Layout, trans, n, k);

    // C_1 = A_1 * A_0 // grade school
    // std::cerr << "\tmultiplying a[" << a_start + len2 << ":" <<
    // a_start+2*len2 << "] * a[" << a_start << ":" << a_start+len2 << "] -> c["
    // << c_start + (3*clen2+2)<< ":" <<  c_start + (3*clen2+2)+2*clen2 + 1 <<
    // "]\n";
    gradeschool_cpu(a_start + len2, a_start, c_start + (3 * clen2 + 2), c_max,
                    Layout, trans == CblasTrans ? CblasNoTrans : CblasTrans,
                    trans, n, n, k);

    // move over C_1, deal with the overlap
    // std::cerr << "\tadding c[" << c_start+(3*clen2+2) << ":" <<
    // c_start+(3*clen2+2)+2*clen2+1 << "] -> c[" <<(c_start+clen2+1) << ":" <<
    // (c_start+3*clen2+2) <<"]\n";
    cblas_daxpy(n * n * (2 * clen2 + 1), 1.0,
                c_double_array + n * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + n * n * (c_start + clen2 + 1), 1);

    // C_1 += (A_1 * A_0)T = A_0 * A_1
    // std::cerr << "TRANSPOSING\n";
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
    // std::cerr << "TRANSPOSED\n";
    // std::cerr << "\tadding c[" << c_start+(3*clen2+2) << ":" <<
    // c_start+(3*clen2+2)+2*clen2+1 << "]T -> c[" <<(c_start+clen2+1) << ":" <<
    // (c_start+3*clen2+2) <<"]\n";
    cblas_daxpy(n * n * (2 * clen2 + 1), 1.0,
                c_double_array + n * n * (c_start + (3 * clen2 + 2)), 1,
                c_double_array + n * n * (c_start + clen2 + 1), 1);

    // clean up
    memset(c_double_array + n * n * (c_start + (3 * clen2 + 2)), 0,
           n * n * (2 * clen2 + 1) * sizeof(double));
  }
}
// syrk gradeschool
// IMPORTANT STEP: please make sure that a_double_array = b_double_array ahead
// of time
void mpmat::gradeschool_gpu(const int &a_start, const int &c_start,
                            const int &c_max, CBLAS_ORDER Layout,
                            CBLAS_TRANSPOSE trans, const int &n, const int &k,
                            int device, const double alpha, const double beta) {
  cudaSetDevice(device);
  int start = 2 * a_start;
  int diff = c_max - start;
  // std::cerr << "sgradeschool (n,k)=(" << n << "," << k << ") " << start << "
  // " << c_start << " " << c_max << "\n";
  if (diff < 2) { // base case, just multiply
    cublasDsyrk(handles[device], CUBLAS_FILL_MODE_UPPER,
                (Layout == CblasRowMajor) != (trans == CblasTrans)
                    ? CUBLAS_OP_N
                    : CUBLAS_OP_T,
                n, k, &alpha, d_a[device] + k * n * a_start,
                (Layout == CblasRowMajor) != (trans == CblasTrans) ? n : k,
                &beta, d_c[device] + n * n * c_start, n);
  }
  // REVIEW
  else { // if we don't need all four, then just do grade school
    int len2 = pow(2, ceil(log2(diff)) - 1) + .1;
    int clen2 =
        len2 -
        1; // the "fundamental length" of c one level below (post squishing)

    // C_0 = A_1 * A_0 // grade school
    gradeschool_gpu(a_start, a_start + len2, c_start, c_max, Layout,
                    trans == CblasTrans ? CblasNoTrans : CblasTrans, trans, n,
                    n, k);

    // C_1 += (A_1 * A_0)T = A_0 * A_1
    for (int i = 0; i < 2 * clen2 + 1; ++i)
      cublasDgeam(handles[device], CUBLAS_OP_T, CUBLAS_OP_N, n, n, &alpha,
                  d_c[device] + n * n * (c_start + i), n, &beta,
                  d_c[device] + n * n * (c_start + (3 * clen2 + 2) + i), n,
                  d_c[device] + n * n * (c_start + (3 * clen2 + 2) + i), n);

    // move over C_0, deal with the overlap
    cublasDaxpy(handles[device], n * n * (2 * clen2 + 1), &alpha,
                d_c[device] + n * n * (c_start), 1,
                d_c[device] + n * n * (c_start + 3 * clen2 + 2), 1);

    cudaMemset(d_c[device] + n * n * (c_start), 0,
               n * n * (2 * clen2 + 1) * sizeof(double));

    // C_0 = A_0 * A_0 // karatsuba
    karatsuba_gpu(a_start, c_start, c_max, Layout, trans, n, k);

    // move over C_1, deal with the overlap
    cublasDaxpy(handles[device], n * n * (2 * clen2 + 1), &alpha,
                d_c[device] + n * n * (c_start + (3 * clen2 + 2)), 1,
                d_c[device] + n * n * (c_start + clen2 + 1), 1);

    // clean up
    cudaMemset(d_c[device] + n * n * (c_start + (3 * clen2 + 2)), 0,
               n * n * (2 * clen2 + 1) * sizeof(double));
  }
}
#endif
// given a length of a trinary tree, condense it into the nearly-binary tree
// that we need specifically a 3^n tree into a 2^(n+1)-1
void mpmat::treecondense(double *c, int size, int l) {
  if (l < 9) {
    return;
  } else {
    int len3 = l / 9; // len3 = 3^(n-2), the length of the leaves two levels
                      // down. the len3 at the next level is 3*len3
    int n = log(l) / log(3) +
            .1; // tree level // the epsilon is because floating point sucks
    int len2 =
        (1 << (n - 1)) - 1; // the length of leaves two levels down post-squish
                            // the len2 at the next level is 2*len2 + 1

    // std::cout << "given l=" << l << ": len3=" << len3 << ", n=" << n << ",
    // len2=" << len2 << "\n";

    // recurse! but carefully so as to not mess up the pointers
    treecondense(c + size * len3 * 6, size, l / 3);
    treecondense(c + size * len3 * 3, size, l / 3);
    treecondense(c, size, l / 3);

    // if (l < 243){

    // now, a delicate operation: stitching together the newly squished answers.
    // Assuming that the previous level properly squished each answer,
    // the size went from 3^(n-1) to 2^n - 1.
    // C_02 += C_10
    // std::cout << "\t\t\tadd array of length " << len2 << " from " << len3*3
    // << " to " << (len2+1) << "\n";
    cblas_daxpy(size * len2, 1.0, c + size * len3 * 3, 1, c + size * (len2 + 1),
                1);
    // C_12 += C_20
    // std::cout << "\t\t\tadd array of length " << len2 << " from " << len3*6
    // << " to " << (len3*3+len2+1) << "\n";
    cblas_daxpy(size * len2, 1.0, c + size * len3 * 6, 1,
                c + size * (len3 * 3 + len2 + 1), 1);

    // Fill in the gap between C_0 and C_1 without overlaps
    // i.e. overwrite the gap and where C_10 was
    // std::cout << "\t\t\tcopy array of length " << len2 << " from " <<
    // (len3*3+len2) << " to " << (2*len2+1) << "\n";
    cblas_dcopy(size * len2, c + size * (len3 * 3 + len2), 1,
                c + size * (2 * len2 + 1), 1);
    // Fill in the rest
    // std::cout << "\t\t\tcopy array of length " << 1 << " from " <<
    // (len3*3+2*len2) << " to " << (3*len2+1) << "\n";
    cblas_dcopy(size, c + size * (len3 * 3 + 2 * len2), 1,
                c + size * (3 * len2 + 1), 1);

    // get the last chunk. 3*len2+2 + len2+1 == 2^(n+1) - 1
    // std::cout << "\t\t\tcopy array of length " << len2+1 << " from " <<
    // len3*6+len2 << " to " << 3*len2+2 << "\n";
    cblas_dcopy(size * (len2 + 1), c + size * (len3 * 6 + len2), 1,
                c + size * (3 * len2 + 2), 1);

    // clean up a bit, juust in case
    // memset((c+size*((int)pow(2,n+1)-1)), 0, (l*size -pow(2,n+1)));
    //}
    return;
  }
}

void mpmat::karatsuba_bc(const int &a_start, const int &b_start,
                         const int &c_start, const int &c_max,
                         CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                         CBLAS_TRANSPOSE transb, const int &m, const int &n,
                         const int &k, const double alpha, const double beta) {

  int start = a_start + b_start;
  int diff = c_max - start;
  std::cout << "\tkaratsuba_bc start=(" << a_start << ", " << b_start
            << ") l=" << c_max << "\n";
  if (diff <= 2) { // base case, just multiply
    cblas_dgemm(
        Layout, transa, transb, m, n, k, alpha,
        a_double_array + k * m * a_start, Layout == CblasRowMajor ? k : m,
        b_double_array + k * n * b_start, Layout == CblasRowMajor ? n : k, beta,
        c_double_array + m * n * c_start, Layout == CblasRowMajor ? n : m);
    // std::cout << "\t\tc = " << c_double_array[m*n*c_start] << "\n";
  } else { // if we need all four quadrants, do full karatsuba
    int len2 = pow(2, ceil(log2(diff)) - 2);
    int len3 = pow(3, ceil(log2(diff)) - 2);
    // C_0 = A_0 * B_0
    karatsuba_bc(a_start, b_start, c_start, 2 * len2 + start, Layout, transa,
                 transb, m, n, k);
    // C_2 = A_1 * B_1 // conditional karatsuba
    if (diff <=
        3 * len2) // if we need half or less of C_2, then just use grade school
      gradeschool_bc(a_start + len2, b_start + len2, c_start + 2 * len3, c_max,
                     Layout, transa, transb, m, n, k);
    else // otherwise, we need all quadrants and therefore karatsuba
      karatsuba_bc(a_start + len2, b_start + len2, c_start + 2 * len3, c_max,
                   Layout, transa, transb, m, n, k);
    // C_1 = A_1 * B_0 // full karatsuba
    karatsuba_bc(a_start + len2, b_start, c_start + len3, 2 * len2 + start,
                 Layout, transa, transb, m, n, k);
    // C_1 = A_0 * B_1 // full karatsuba
    karatsuba_bc(a_start, b_start + len2, c_start + len3, 2 * len2 + start,
                 Layout, transa, transb, m, n, k);

    // std::cout << "\n";
  }
}

void mpmat::gradeschool_bc(const int &a_start, const int &b_start,
                           const int &c_start, const int &c_max,
                           CBLAS_ORDER Layout, CBLAS_TRANSPOSE transa,
                           CBLAS_TRANSPOSE transb, const int &m, const int &n,
                           const int &k, const double alpha,
                           const double beta) {
  int start = a_start + b_start;
  int diff = c_max - start;
  // std::cout << "\tgradeschool_bc start=(" << a_start << ", " << b_start << ")
  // l=" << c_max << "\n";
  if (diff < 2) { // base case, just multiply
    // std::cout << "\t\tbasecase'd out\n";
    cblas_dgemm(
        Layout, transa, transb, m, n, k, alpha,
        a_double_array + k * m * a_start, Layout == CblasRowMajor ? k : m,
        b_double_array + k * n * b_start, Layout == CblasRowMajor ? n : k, beta,
        c_double_array + m * n * c_start, Layout == CblasRowMajor ? n : m);
    // std::cout << "\t\tc = " << c_double_array[m*n*c_start] << "\n";
  } else { // if we don't need all four, then just do grade school
    int len2 = pow(2, ceil(log2(diff)) - 1);
    int len3 = pow(3, ceil(log2(diff)) - 1);
    // C_0 = A_0 * B_0 // karatsuba
    karatsuba_bc(a_start, b_start, c_start, c_max, Layout, transa, transb, m, n,
                 k);
    // C_1 = A_1 * B_0 // grade school
    gradeschool_bc(a_start + len2, b_start, c_start + len3, c_max, Layout,
                   transa, transb, m, n, k);
    // C_2 = A_0 * B_1 // grade school // stored in C_2 temporarily!
    gradeschool_bc(a_start, b_start + len2, c_start + 2 * len3, c_max, Layout,
                   transa, transb, m, n, k);
    // std::cout << "\n";
    // C_1 += C_2
    cblas_daxpy(m * n * len3, 1.0,
                c_double_array + m * n * (c_start + 2 * len3), 1,
                c_double_array + m * n * (c_start + len3), 1);
    // C_2 = 0
    memset(c_double_array + m * n * (c_start + 2 * len3), 0,
           m * n * len3 * sizeof(double));
  }
}
