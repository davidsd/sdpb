//
// Created by Rajeev Erramilli on 12/1/17
//

#include "../Timers.h"
#include "mpmat.h"
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>

#include <bitset>
#include <iostream>

template <typename T> inline T ceil_div(T a, T b) {
  return a / b + ((a % b != 0) & (a > 0));
}

double *randomDoubleVector(int n) {
  auto out = new double[n];
  std::uniform_real_distribution<double> unif(1, 5);
  std::random_device rd;
#pragma omp parallel for
  for (int i = 0; i < n; ++i)
    out[i] = floor(unif(rd));
  return out;
}

bool compareSymmMatrices(double *x, double *y, int n, int l, bool lower,
                         bool rowmaj) {
  if (lower != !rowmaj) {
    for (int i = 0; i < l; ++i) {
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k <= j; ++k) {
          if (x[n * n * i + n * j + k] != y[n * n * i + n * j + k]) {
            std::cerr << "Error: matrices disagree at (" << i << "," << j << ","
                      << k << ").\n";
            return false;
          }
        }
      }
    }
  } else {
    for (int i = 0; i < l; ++i) {
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k <= j; ++k) {
          if (x[n * n * i + n * (n - j - 1) + (n - k)] !=
              y[n * n * i + n * (n - j - 1) + (n - k)]) {
            std::cerr << "Error: matrices disagree at (" << i << "," << j << ","
                      << k << ").\n";
            return false;
          }
        }
      }
    }
  }
  return true;
}

void mpmat::mpmat_conversion_test(int i, int f, int d) {
  mpf_class in = i;
  mpf_set_default_prec(mpf_get_default_prec() + 256);
  mpf_class out;
  mpf_set_default_prec(mpf_get_default_prec() - 256);
  int mpmat_limb = (MPMAT_DOUBLE_MANT_IMPLICIT - 13) / 2;
  int mpmat_size = ceil_div(
      abs(in.get_mpf_t()->_mp_prec + 1) * mp_bits_per_limb, mpmat_limb);
  int exp = in.get_mpf_t()->_mp_exp * mp_bits_per_limb;
  double *arr = new double[mpmat_size];
  for (; in < f; in *= d) {
    std::cout << "Testing conversion of " << in << "...\n";
    mpmatConvertGMPToDouble(in, arr, mpmat_size, mpmat_limb, exp);
    mpmatConvertDoubleToGMP(out, arr, mpmat_size, mpmat_limb, exp);
    compare_mpf_bits(in, out);
    memset(arr, 0, mpmat_size * sizeof(double));
    std::cout << "\n\n";
  }
}

bool mpmat::karatsuba_test(int m, int n, int k, int l) {
  double *tmp_a = a_double_array, *tmp_b = b_double_array,
         *tmp_c = c_double_array;
  a_double_array = randomDoubleVector(m * k * l);
  double *a2_double_array = new double[m * k * l];
  b_double_array = randomDoubleVector(n * k * l);
  double *b2_double_array = new double[n * k * l];
  c_double_array = new double[m * n * (6 * l - (int)log2(l) - 2)];
  double *c2_double_array = new double[m * n * l];

  std::copy(a_double_array, a_double_array + m * k * l, a2_double_array);
  std::copy(b_double_array, b_double_array + n * k * l, b2_double_array);

  // std::cout << "\ta_double_array: ";
  // for (int i = 0; i < m*k*l; ++i)
  //   std::cout << a_double_array[i] << " ";
  // std::cout << "\n";

  // std::cout << "\tb_double_array: ";
  // for (int i = 0; i < n*k*l; ++i)
  //   std::cout << b_double_array[i] << " ";
  // std::cout << "\n";
  memset(c2_double_array, 0, m * n * l * sizeof(double));
  memset(c_double_array, 0,
         m * n * (6 * l - (int)log2(l) - 2) * sizeof(double));

  timers[std::string("gradeschool.l=") + std::to_string(l)].start();
  for (int i = 0; i < l; i++)
    for (int j = 0; j <= i; j++)
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1,
                  a_double_array + k * m * j, k,
                  b_double_array + (i - j) * k * n, n, 1,
                  c2_double_array + i * m * n, n);
  timers[std::string("gradeschool.l=") + std::to_string(l)].stop();
  // std::cout << a_double_array[0] << " times " << b_double_array[0] << "
  // equals " << c2_double_array[0] << ".\n";

  // gradeschool_bc(0, 0, 0, l, CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n,
  // k);
  // // std::cout << "\tbefore squish: ";
  // // for (int i = 0; i < m*n*(int)pow(3,ceil(log2(l))); ++i)
  // //   std::cout << c_double_array[i] << " ";
  // // std::cout << "\n";

  // treecondense(c_double_array, m*n, (int)pow(3,ceil(log2(l))));
  // // std::cout << "\tafter squish: ";
  // // for (int i = 0; i < m*n*l; ++i)
  // //     std::cout << c_double_array[i] << " ";
  // // std::cout << "\n";

  memset(c_double_array, 0,
         m * n * (6 * l - (int)log2(l) - 2) * sizeof(double));
  timers[std::string("karatsuba.l=") + std::to_string(l)].start();
  karatsuba_generic(l, CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k);

  // std::cout << "\tbefore squish: ";
  //  for (int i = 0; i < m*n*(int)pow(3,ceil(log2(l))); ++i)
  //    std::cout << c_double_array[i] << " ";
  // std::cout << "\n";

  // std::cout << "condensing a tree with initial length: " <<
  // (int)pow(3,ceil(log2(l))) << "\n";

  // treecondense(c_double_array, m*n, (int)pow(3,ceil(log2(l))));
  timers[std::string("karatsuba.l=") + std::to_string(l)].stop();

  // std::cout << "\tafter squish: ";
  // for (int i = 0; i < m*n*l; ++i)
  //     std::cout << c_double_array[i] << " ";
  // std::cout << "\n";

  bool result = true;

  if (!std::equal(a_double_array, a_double_array + m * k * l,
                  a2_double_array)) {
    std::cerr << "a_double_array was modified\n";
    result = false;
    for (int i = 0; i < m * k * l; ++i)
      std::cout << a_double_array[i] << " ";
    std::cout << "\n";
    for (int i = 0; i < m * k * l; ++i)
      std::cout << a2_double_array[i] << " ";
    std::cout << "\n\n";
  }
  if (!std::equal(b_double_array, b_double_array + n * k * l,
                  b2_double_array)) {
    std::cerr << "b_double_array was modified\n";
    result = false;
    for (int i = 0; i < n * k * l; ++i)
      std::cout << b_double_array[i] << " ";
    std::cout << "\n";
    for (int i = 0; i < n * k * l; ++i)
      std::cout << b2_double_array[i] << " ";
    std::cout << "\n\n";
  }

  if (!std::equal(c2_double_array, c2_double_array + n * m * l,
                  c_double_array)) {
    std::cerr << "c_double_array is incorrect\n";
    result = false;
    for (int i = 0; i < m * n * (6 * l - (int)log2(l) - 2); ++i)
      std::cout << c_double_array[i] << " ";
    std::cout << "\n";
    for (int i = 0; i < m * n * l; ++i)
      std::cout << c2_double_array[i] << " ";
    std::cout << "\n";
    // print_mpmat_double_array(c_double_array,m*n*(int)pow(3,ceil(log2(l))));
    // print_mpmat_double_array(c2_double_array,m*n*l);
  }

  delete[] a_double_array;
  delete[] b_double_array;
  delete[] c_double_array;
  delete[] a2_double_array;
  delete[] b2_double_array;
  delete[] c2_double_array;

  a_double_array = tmp_a;
  b_double_array = tmp_b;
  c_double_array = tmp_c;

  return result;
}

bool mpmat::symm_karatsuba_test(int n, int k, int l) {
  timers[std::string("setup.l=") + std::to_string(n)].start();
  double *tmp_a = a_double_array, *tmp_c = c_double_array;
  a_double_array = randomDoubleVector(n * k * l);
  double *a2_double_array = new double[n * k * l];
  c_double_array = new double[n * n * (6 * l + 2)];
  double *b2_double_array = new double[n * n * l];
  double *c2_double_array = new double[n * n * l];

  std::copy(a_double_array, a_double_array + n * k * l, a2_double_array);

  // std::cout << "\ta_double_array: ";
  // for (int i = 0; i < n*k*l; ++i)
  //   std::cout << a_double_array[i] << " ";
  // std::cout << "\n";

  // std::cout << "\tb_double_array: ";
  // for (int i = 0; i < n*k*l; ++i)
  //   std::cout << b_double_array[i] << " ";
  // std::cout << "\n";
  memset(c2_double_array, 0, n * n * l * sizeof(double));
  memset(b2_double_array, 0, n * n * l * sizeof(double));
  memset(c_double_array, 0,
         n * n * (3 * l - (int)log2(l) - 1) * sizeof(double));
  timers[std::string("setup.l=") + std::to_string(n)].stop();

  timers[std::string("gradeschool.l=") + std::to_string(n)].start();
  for (int i = 0; i < l; i++) {
    for (int j = 0; j < i / 2 + i % 2; j++) {
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, k, 1.0,
                  a_double_array + k * n * j, k,
                  a_double_array + (i - j) * k * n, k, 1.0,
                  c2_double_array + i * n * n, n);
    }
#ifdef HAVE_MKL_H
    mkl_domatcopy('c', 't',
#else
    cblas_domatcopy(CblasColMajor, CblasTrans,
#endif
                  n, n, 1, c2_double_array + i * n * n, n,
                  b2_double_array + i * n * n, n);
    cblas_daxpy(n * n, 1.0, b2_double_array + i * n * n, 1,
                c2_double_array + i * n * n, 1);
    // if significance of result is even, calculate the symmetric part
    if (i % 2 == 0)
      cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, n, k, 1.0,
                  a_double_array + k * n * (i / 2), k, 1.0,
                  c2_double_array + i * n * n, n);
  }
  timers[std::string("gradeschool.l=") + std::to_string(n)].stop();
  // std::cout << a_double_array[0] << " times " << b_double_array[0] << "
  // equals " << c2_double_array[0] << ".\n";

  // gradeschool_bc(0, 0, 0, l, CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n,
  // k);
  // // std::cout << "\tbefore squish: ";
  // // for (int i = 0; i < m*n*(int)pow(3,ceil(log2(l))); ++i)
  // //   std::cout << c_double_array[i] << " ";
  // // std::cout << "\n";

  // treecondense(c_double_array, m*n, (int)pow(3,ceil(log2(l))));
  // // std::cout << "\tafter squish: ";
  // // for (int i = 0; i < m*n*l; ++i)
  // //     std::cout << c_double_array[i] << " ";
  // // std::cout << "\n";

  memset(c_double_array, 0,
         n * n * (6 * l - (int)log2(l) - 1) * sizeof(double));
  timers[std::string("karatsuba.l=") + std::to_string(n)].start();
  karatsuba_symmetric(l, CblasRowMajor, CblasTrans, n, k);
  // std::cerr << "done with karatsuba test\n";

  // std::cout << "\tbefore squish: ";
  //  for (int i = 0; i < m*n*(int)pow(3,ceil(log2(l))); ++i)
  //    std::cout << c_double_array[i] << " ";
  // std::cout << "\n";

  // std::cout << "condensing a tree with initial length: " <<
  // (int)pow(3,ceil(log2(l))) << "\n";

  // treecondense(c_double_array, m*n, (int)pow(3,ceil(log2(l))));
  timers[std::string("karatsuba.l=") + std::to_string(n)].stop();

  // std::cout << "\tafter squish: ";
  // for (int i = 0; i < m*n*l; ++i)
  //     std::cout << c_double_array[i] << " ";
  // std::cout << "\n";

  bool result = true;

  if (!std::equal(a_double_array, a_double_array + n * k * l,
                  a2_double_array)) {
    std::cerr << "a_double_array was modified\n";
    result = false;
    for (int i = 0; i < n * k * l; ++i)
      std::cout << a_double_array[i] << " ";
    std::cout << "\n";
    for (int i = 0; i < n * k * l; ++i)
      std::cout << a2_double_array[i] << " ";
    std::cout << "\n\n";
  }
  if (!compareSymmMatrices(c_double_array, c2_double_array, n, l)) {
    std::cerr << "c_double_array is incorrect\n";
    result = false;
    for (int i = 0; i < n * n * (3 * l - (int)log2(l) - 1); ++i)
      std::cout << c_double_array[i]
                << (i % (n * n) == n * n - 1 ? "|"
                                             : (i % n == n - 1 ? "," : " "));
    std::cout << "\n";
    for (int i = 0; i < n * n * l; ++i)
      std::cout << c2_double_array[i]
                << (i % (n * n) == n * n - 1 ? "|"
                                             : (i % n == n - 1 ? "," : " "));
    std::cout << "\n";
  }

  timers[std::string("cleanup.l=") + std::to_string(n)].start();

  delete[] a_double_array;
  delete[] c_double_array;
  delete[] a2_double_array;
  delete[] b2_double_array;
  delete[] c2_double_array;

  a_double_array = tmp_a;
  c_double_array = tmp_c;

  timers[std::string("cleanup.l=") + std::to_string(n)].stop();
  return result;
}
#ifdef __SDPB_CUDA__
bool mpmat::karatsuba_test_gpu(int m, int n, int k, int l) {
  double *tmp_a = a_double_array, *tmp_b = b_double_array,
         *tmp_c = c_double_array;
  // double * tmp_d_a = d_a[0], * tmp_d_b = d_b[0], * tmp_d_c = d_c[0];
  a_double_array = randomDoubleVector(m * k * l);
  // cudaMalloc(d_a,m*k*l*sizeof(mpmat_double));
  double *a2_double_array = new double[m * k * l];
  b_double_array = randomDoubleVector(n * k * l);
  // cudaMalloc(d_b,n*k*l*sizeof(mpmat_double));
  double *b2_double_array = new double[n * k * l];
  c_double_array = new double[m * n * (6 * l - (int)log2(l) - 2)];
  // cudaMalloc(d_c,m*n*(6*l - (int)log2(l) - 2)*sizeof(mpmat_double));
  double *c2_double_array = new double[m * n * l];

  std::copy(a_double_array, a_double_array + m * k * l, a2_double_array);
  std::copy(b_double_array, b_double_array + n * k * l, b2_double_array);

  std::cout << "\ta_double_array: ";
  for (int i = 0; i < m * k * l; ++i)
    std::cout << a_double_array[i] << " ";
  std::cout << "\n";

  std::cout << "\tb_double_array: ";
  for (int i = 0; i < n * k * l; ++i)
    std::cout << b_double_array[i] << " ";
  std::cout << "\n";
  memset(c2_double_array, 0, m * n * l * sizeof(double));
  memset(c_double_array, 0,
         m * n * (6 * l - (int)log2(l) - 2) * sizeof(double));

  timers[std::string("gradeschool.l=") + std::to_string(l)].start();
  for (int i = 0; i < l; i++)
    for (int j = 0; j <= i; j++)
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1,
                  a_double_array + k * m * j, k,
                  b_double_array + (i - j) * k * n, n, 1,
                  c2_double_array + i * m * n, n);
  timers[std::string("gradeschool.l=") + std::to_string(l)].stop();
  // std::cout << a_double_array[0] << " times " << b_double_array[0] << "
  // equals " << c2_double_array[0] << ".\n";

  // gradeschool_bc(0, 0, 0, l, CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n,
  // k);
  // // std::cout << "\tbefore squish: ";
  // // for (int i = 0; i < m*n*(int)pow(3,ceil(log2(l))); ++i)
  // //   std::cout << c_double_array[i] << " ";
  // // std::cout << "\n";

  // treecondense(c_double_array, m*n, (int)pow(3,ceil(log2(l))));
  // // std::cout << "\tafter squish: ";
  // // for (int i = 0; i < m*n*l; ++i)
  // //     std::cout << c_double_array[i] << " ";
  // // std::cout << "\n";

  // memset(c_double_array, 0, m*n*(6*l - (int)log2(l) - 2)*sizeof(double));

  timers[std::string("karatsuba.l=") + std::to_string(l)].start();
  // cudaMemcpy(d_a[0],a_double_array,m*k*l*sizeof(mpmat_double),cudaMemcpyHostToDevice);
  // cudaMemcpy(d_b[0],b_double_array,n*k*l*sizeof(mpmat_double),cudaMemcpyHostToDevice);

  // cudaMemcpy(a_double_array,d_a[0],m*k*l*sizeof(mpmat_double),cudaMemcpyDeviceToHost);
  // cudaMemcpy(b_double_array,d_b[0],n*k*l*sizeof(mpmat_double),cudaMemcpyDeviceToHost);
  //   std::cout << "\ta_double_array: ";
  // for (int i = 0; i < m*k*l; ++i)
  //   std::cout << a_double_array[i] << " ";
  // std::cout << "\n";

  // std::cout << "\tb_double_array: ";
  // for (int i = 0; i < n*k*l; ++i)
  //   std::cout << b_double_array[i] << " ";
  // std::cout << "\n";

  // cudaMemset(d_c[0],0,m*n*(6*l - (int)log2(l) - 2)*sizeof(double));

  karatsuba_generic(l, CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, true);

  // cudaMemcpy(c_double_array,d_c[0],m*n*(6*l - (int)log2(l) -
  // 2)*sizeof(mpmat_double),cudaMemcpyDeviceToHost);

  // std::cout << "\tbefore squish: ";
  //  for (int i = 0; i < m*n*(int)pow(3,ceil(log2(l))); ++i)
  //    std::cout << c_double_array[i] << " ";
  // std::cout << "\n";

  // std::cout << "condensing a tree with initial length: " <<
  // (int)pow(3,ceil(log2(l))) << "\n";

  // treecondense(c_double_array, m*n, (int)pow(3,ceil(log2(l))));
  timers[std::string("karatsuba.l=") + std::to_string(l)].stop();

  // std::cout << "\tafter squish: ";
  // for (int i = 0; i < m*n*l; ++i)
  //     std::cout << c_double_array[i] << " ";
  // std::cout << "\n";

  bool result = true;

  if (!std::equal(a_double_array, a_double_array + m * k * l,
                  a2_double_array)) {
    std::cerr << "a_double_array was modified\n";
    result = false;
    for (int i = 0; i < m * k * l; ++i)
      std::cout << a_double_array[i] << " ";
    std::cout << "\n";
    for (int i = 0; i < m * k * l; ++i)
      std::cout << a2_double_array[i] << " ";
    std::cout << "\n\n";
  }
  if (!std::equal(b_double_array, b_double_array + n * k * l,
                  b2_double_array)) {
    std::cerr << "b_double_array was modified\n";
    result = false;
    for (int i = 0; i < n * k * l; ++i)
      std::cout << b_double_array[i] << " ";
    std::cout << "\n";
    for (int i = 0; i < n * k * l; ++i)
      std::cout << b2_double_array[i] << " ";
    std::cout << "\n\n";
  }

  if (!std::equal(c2_double_array, c2_double_array + n * m * l,
                  c_double_array)) {
    std::cerr << "c_double_array is incorrect\n";
    result = false;
    for (int i = 0; i < m * n * (6 * l - (int)log2(l) - 2); ++i)
      std::cout << c_double_array[i] << " ";
    std::cout << "\n";
    for (int i = 0; i < m * n * l; ++i)
      std::cout << c2_double_array[i] << " ";
    std::cout << "\n";
    // print_mpmat_double_array(c_double_array,m*n*(int)pow(3,ceil(log2(l))));
    // print_mpmat_double_array(c2_double_array,m*n*l);
  }

  delete[] a_double_array;
  delete[] b_double_array;
  delete[] c_double_array;
  delete[] a2_double_array;
  delete[] b2_double_array;
  delete[] c2_double_array;

  // cudaFree(d_a[0]);
  // cudaFree(d_b[0]);
  // cudaFree(d_c[0]);

  a_double_array = tmp_a;
  b_double_array = tmp_b;
  c_double_array = tmp_c;
  // d_a[0] = tmp_d_a;
  // d_b[0] = tmp_d_b;
  // d_c[0] = tmp_d_c;

  return result;
}

bool mpmat::symm_karatsuba_test_gpu(int n, int k, int l) {
  timers[std::string("setup.l=") + std::to_string(n)].start();
  double *tmp_a = a_double_array, *tmp_c = c_double_array;
  // double * tmp_d_a = d_a[0], * tmp_d_b = d_b[0], * tmp_d_c = d_c[0];
  a_double_array = randomDoubleVector(n * k * l);
  // cudaMalloc(d_a,n*k*l*sizeof(mpmat_double));
  // cudaMalloc(d_b,n*k*l*sizeof(mpmat_double));
  double *a2_double_array = new double[n * k * l];
  double *b2_double_array = new double[n * n * l];
  c_double_array = new double[n * n * (6 * l - (int)log2(l) - 2)];
  // cudaMalloc(d_c,n*n*(6*l - (int)log2(l) - 2)*sizeof(mpmat_double));
  double *c2_double_array = new double[n * n * l];

  std::copy(a_double_array, a_double_array + n * k * l, a2_double_array);
  // std::copy(b_double_array,b_double_array+n*k*l,b2_double_array);

  std::cout << "\ta_double_array: ";
  for (int i = 0; i < n * k * l; ++i)
    std::cout << a_double_array[i] << " ";
  std::cout << "\n";

  // std::cout << "\tb_double_array: ";
  // for (int i = 0; i < n*k*l; ++i)
  //   std::cout << b_double_array[i] << " ";
  // std::cout << "\n";
  memset(c2_double_array, 0, n * n * l * sizeof(double));
  memset(c_double_array, 0,
         n * n * (3 * l - (int)log2(l) - 2) * sizeof(double));
  timers[std::string("setup.l=") + std::to_string(n)].stop();

  std::cerr << "starting gradeschool test\n";
  timers[std::string("gradeschool.l=") + std::to_string(n)].start();
  for (int i = 0; i < l; i++) {
    for (int j = 0; j < i / 2 + i % 2; j++) {
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, k, 1.0,
                  a_double_array + k * n * j, k,
                  a_double_array + (i - j) * k * n, k, 1.0,
                  c2_double_array + i * n * n, n);
    }
#ifdef HAVE_MKL_H
    mkl_domatcopy('c', 't',
#else
    cblas_domatcopy(CblasColMajor, CblasTrans,
#endif
                  n, n, 1, c2_double_array + i * n * n, n,
                  b2_double_array + i * n * n, n);
    cblas_daxpy(n * n, 1.0, b2_double_array + i * n * n, 1,
                c2_double_array + i * n * n, 1);
    // if significance of result is even, calculate the symmetric part
    if (i % 2 == 0)
      cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, n, k, 1.0,
                  a_double_array + k * n * (i / 2), k, 1.0,
                  c2_double_array + i * n * n, n);
  }
  timers[std::string("gradeschool.l=") + std::to_string(n)].stop();
  // std::cout << a_double_array[0] << " times " << b_double_array[0] << "
  // equals " << c2_double_array[0] << ".\n";

  // gradeschool_bc(0, 0, 0, l, CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n,
  // k);
  // // std::cout << "\tbefore squish: ";
  // // for (int i = 0; i < m*n*(int)pow(3,ceil(log2(l))); ++i)
  // //   std::cout << c_double_array[i] << " ";
  // // std::cout << "\n";

  // treecondense(c_double_array, m*n, (int)pow(3,ceil(log2(l))));
  // // std::cout << "\tafter squish: ";
  // // for (int i = 0; i < m*n*l; ++i)
  // //     std::cout << c_double_array[i] << " ";
  // // std::cout << "\n";

  std::cerr << "starting karatsuba test\n";
  memset(c_double_array, 0,
         n * n * (6 * l - (int)log2(l) - 2) * sizeof(double));
  timers[std::string("karatsuba.l=") + std::to_string(n)].start();

  // cudaMemcpy(d_a[0],a_double_array,n*k*l*sizeof(mpmat_double),cudaMemcpyHostToDevice);
  // std::cout << "\ta_double_array: ";
  // for (int i = 0; i < n*k*l; ++i)
  //   std::cout << a_double_array[i] << " ";
  // std::cout << "\n";

  // cudaMemset(d_c[0],0,n*n*(6*l - (int)log2(l) - 2)*sizeof(double));

  karatsuba_symmetric(l, CblasRowMajor, CblasTrans, n, k, true);
  // cudaMemcpy(c_double_array,d_c[0],n*n*(6*l - (int)log2(l) -
  // 2)*sizeof(mpmat_double),cudaMemcpyDeviceToHost);

  // std::cerr << "done with karatsuba test\n";

  // std::cout << "\tbefore squish: ";
  //  for (int i = 0; i < m*n*(int)pow(3,ceil(log2(l))); ++i)
  //    std::cout << c_double_array[i] << " ";
  // std::cout << "\n";

  // std::cout << "condensing a tree with initial length: " <<
  // (int)pow(3,ceil(log2(l))) << "\n";

  // treecondense(c_double_array, m*n, (int)pow(3,ceil(log2(l))));
  timers[std::string("karatsuba.l=") + std::to_string(n)].stop();

  // std::cout << "\tafter squish: ";
  // for (int i = 0; i < m*n*l; ++i)
  //     std::cout << c_double_array[i] << " ";
  // std::cout << "\n";

  bool result = true;

  if (!std::equal(a_double_array, a_double_array + n * k * l,
                  a2_double_array)) {
    std::cerr << "a_double_array was modified\n";
    result = false;
    for (int i = 0; i < n * k * l; ++i)
      std::cout << a_double_array[i] << " ";
    std::cout << "\n";
    for (int i = 0; i < n * k * l; ++i)
      std::cout << a2_double_array[i] << " ";
    std::cout << "\n\n";
  }
  if (!compareSymmMatrices(
          c_double_array, c2_double_array, n,
          l)) { //(!std::equal(c2_double_array,c2_double_array+n*n*l,c_double_array)){
    std::cerr << "c_double_array is incorrect\n";
    result = false;
    for (int i = 0; i < n * n * (3 * l - (int)log2(l) - 1); ++i)
      std::cout << c_double_array[i]
                << (i % (n * n) == n * n - 1 ? "|"
                                             : (i % n == n - 1 ? "," : " "));
    std::cout << "\n";
    for (int i = 0; i < n * n * l; ++i)
      std::cout << c2_double_array[i]
                << (i % (n * n) == n * n - 1 ? "|"
                                             : (i % n == n - 1 ? "," : " "));
    std::cout << "\n";
    // print_mpmat_double_array(c_double_array,m*n*(int)pow(3,ceil(log2(l))));
    // print_mpmat_double_array(c2_double_array,m*n*l);
  }

  timers[std::string("cleanup.l=") + std::to_string(n)].start();

  delete[] a_double_array;
  delete[] c_double_array;
  delete[] a2_double_array;
  delete[] b2_double_array;
  delete[] c2_double_array;

  // cudaFree(d_a[0]);
  // cudaFree(d_b[0]);
  // cudaFree(d_c[0]);

  a_double_array = tmp_a;
  // b_double_array = tmp_b;
  c_double_array = tmp_c;
  // d_a[0] = tmp_d_a;
  // d_b[0] = tmp_d_b;
  // d_c[0] = tmp_d_c;

  timers[std::string("cleanup.l=") + std::to_string(n)].stop();
  return result;
}
#endif
// bool mpmat::base_karatsuba_test(){
//   double * tmp_a = a_double_array, * tmp_b = b_double_array, * tmp_c =
//   c_double_array; a_double_array = new double [1]; a_double_array[0] = 2.0;
//   double * a2_double_array = new double[1];
//   b_double_array = new double [1];
//   b_double_array[0] = 3.0;
//   double * b2_double_array = new double[1];
//   c_double_array = new double [1];
//   double * c2_double_array = new double [1];

//   std::copy(a_double_array,a_double_array+1,a2_double_array);
//   std::cout << a_double_array[0] << " -> " << a2_double_array[0] << "\n";
//   std::copy(b_double_array,b_double_array+1,b2_double_array);
//   std::cout << b_double_array[0] << " -> " << b2_double_array[0] << "\n";

//   int n = 1, k = 1, m = 1, l = 1;

//   for (int i = 0; i < l; i++)
//         for (int j = 0; j <= i; j++)
//             cblas_dgemm(
//                     CblasRowMajor,
//                     CblasNoTrans,
//                     CblasNoTrans,
//                     m,
//                     n,
//                     k,
//                     1,
//                     a_double_array+k*m*j,
//                     k,
//                     b_double_array+(i-j)*k*n,
//                     n,
//                     1,
//                     c2_double_array+i*m*n,
//                     n
//             );
//   std::cout << a_double_array[0] << " times " << b_double_array[0] << "
//   equals " << c2_double_array[0] << ".\n";

//   karatsuba(l, CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k);

//   std::cout << a_double_array[0] << " times " << b_double_array[0] << "
//   equals " << c_double_array[0] << ".\n";

//   treecondense(c_double_array, m*n, l);

//   std::cout << a_double_array[0] << " times " << b_double_array[0] << "
//   equals " << c_double_array[0] << ".\n";

//   bool result = true;

//   if (!std::equal(a_double_array,a_double_array+m*k*l,a2_double_array)){
//     std::cerr << "a_double_array was modified\n";
//     result = false;
//   }
//   if (!std::equal(b_double_array,b_double_array+n*k*l,b2_double_array)){
//     std::cerr << "b_double_array was modified\n";
//     result = false;
//   }

//   if (!std::equal(c2_double_array,c2_double_array+n*k*l,c_double_array)){
//     std::cerr << "c_double_array is incorrect\n";
//     result = false;
//   }

//   delete [] a_double_array;
//   delete [] c_double_array;
//   delete [] a2_double_array;
//   delete [] c2_double_array;

//   a_double_array = tmp_a;
//   c_double_array = tmp_c;

//   return result;
// }
