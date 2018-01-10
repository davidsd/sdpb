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

void mpmat::karatsuba_bc(const int & a_start, const int & b_start, const int & c_start, const int & c_max,
                      CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE transa,
                      CBLAS_TRANSPOSE transb, const int & m, const int & n, const int & k,
                      const double alpha, const double beta){
  
  int start = a_start + b_start;
  int diff = c_max - start;
  //std::cout << "\tkaratsuba_bc start=(" << a_start << ", " << b_start << ") l=" << c_max << "\n";
  if(diff <= 2){ // base case, just multiply
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
    //std::cout << "\t\tc = " << c_double_array[m*n*c_start] << "\n";
  }
  else { //if we need all four quadrants, do full karatsuba
    int len2 = pow(2,ceil(log2(diff))-2);
    int len3 = pow(3,ceil(log2(diff))-2);
    // C_0 = A_0 * B_0 
    karatsuba_bc(a_start, b_start, c_start, 2*len2 + start, Layout, transa, transb, m, n, k);
    // C_2 = A_1 * B_1 // conditional karatsuba
    if (diff <= 3*len2) // if we need half or less of C_2, then just use grade school
      gradeschool_bc(a_start+len2, b_start+len2, c_start+2*len3, c_max, Layout, transa, transb, m, n, k);
    else               // otherwise, we need all quadrants and therefore karatsuba
      karatsuba_bc(a_start+len2, b_start+len2, c_start+2*len3, c_max, Layout, transa, transb, m, n, k);
    // C_1 = A_1 * B_0 // full karatsuba
    karatsuba_bc(a_start+len2, b_start, c_start+len3, 2*len2 + start, Layout, transa, transb, m, n, k);
    // C_1 = A_0 * B_1 // full karatsuba
    karatsuba_bc(a_start, b_start+len2, c_start+len3, 2*len2 + start, Layout, transa, transb, m, n, k);

    //std::cout << "\n";
  }
}

void mpmat::karatsuba(const int & a_start, const int & b_start, const int & c_start, const int & c_max,
                      CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE transa,
                      CBLAS_TRANSPOSE transb, const int & m, const int & n, const int & k,
                      const double alpha, const double beta){
  
  int start = a_start + b_start;
  int diff = c_max - start;
  //std::cout << "\tkaratsuba start=(" << a_start << ", " << b_start << ") l=" << c_max << "\n";
  if(diff <= 2){ // base case, just multiply
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
    //std::cout << "\t\tc = " << c_double_array[m*n*c_start] << "\n";
  }
  else { //if we need all four quadrants, do full karatsuba
    int len2 = pow(2,ceil(log2(diff))-2);
    int len3 = pow(3,ceil(log2(diff))-2);
    // C_0 = A_0 * B_0 // full karatsuba
    karatsuba(a_start, b_start, c_start, 2*len2 + start, Layout, transa, transb, m, n, k);
    // C_2 = A_1 * B_1 // conditional karatsuba
    if (diff <= 3*len2) // if we need half or less of C_2, then just use grade school
      gradeschool(a_start+len2, b_start+len2, c_start+2*len3, c_max, Layout, transa, transb, m, n, k);
    else               // otherwise, we need all quadrants and therefore karatsuba
      karatsuba(a_start+len2, b_start+len2, c_start+2*len3, c_max, Layout, transa, transb, m, n, k);
    // A_0 += A_1
    cblas_daxpy(k*m*len2, 1.0, a_double_array+k*m*(a_start+len2), 1, a_double_array+k*m*a_start, 1);
    // B_0 += B_1
    cblas_daxpy(k*n*len2, 1.0, b_double_array+k*n*(b_start+len2), 1, b_double_array+k*n*b_start, 1);
    // C_1 = A_0 * B_0 // full karatsuba
    karatsuba(a_start, b_start, c_start+len3, 2*len2 + start, Layout, transa, transb, m, n, k);
    //std::cout << "\n";
    // C_1 -= C_0
    cblas_daxpy(m*n*len3, -1.0, c_double_array+m*n*c_start, 1, c_double_array+m*n*(c_start+len3), 1);
    // C_1 -= C_2
    cblas_daxpy(m*n*len3, -1.0, c_double_array+m*n*(c_start+2*len3), 1, c_double_array+m*n*(c_start+len3), 1);
    // A_0 -= A_1 // clean up
    cblas_daxpy(k*m*len2, -1.0, a_double_array+k*m*(a_start+len2), 1, a_double_array+k*m*a_start, 1);
    // B_0 -= B_1 // clean up
    cblas_daxpy(k*n*len2, -1.0, b_double_array+k*n*(b_start+len2), 1, b_double_array+k*n*b_start, 1);
  }
  
}

void mpmat::karatsuba(const int & a_start, const int & c_start, const int & c_max,
                      CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE trans, const int & n, const int & k,
                      const double alpha, const double beta){
  int start = 2*a_start;
  int diff = c_max - start;
  if(diff < 2){ // base case, just multiply
    cblas_dsyrk(CblasColMajor,
                CblasUpper,(Layout == CblasRowMajor) != (trans == CblasTrans) ? CblasNoTrans : CblasTrans,
                n,k,alpha,
                a_double_array+k*n*a_start,
                (Layout == CblasRowMajor) != (trans == CblasTrans) ? n : k,
                beta,
                c_double_array+n*n*c_start,
                n);
  }
  else{ //if we need all four quadrants, do full symmetric karatsuba
    int len2 = pow(2,ceil(log2(diff))-2);
    int len3 = pow(3,ceil(log2(diff))-2);
    // C_0 = A_0 * A_0 // full karatsuba
    karatsuba(a_start, c_start, 2*len2 + start, Layout, trans, n, k);
    // C_2 = A_1 * A_1 // conditional karatsuba
    if (diff <= 3*len2) // if we need half or less of C_2, then just use grade school
      gradeschool(a_start+len2, c_start+2*len3, c_max, Layout, trans, n, k);
    else               // otherwise, we need all quadrants and therefore karatsuba
      karatsuba(a_start+len2, c_start+2*len3, c_max, Layout, trans, n, k);
    // A_0 += A_1
    cblas_daxpy(k*n, 1.0, a_double_array+k*n*(a_start+len2), 1, a_double_array+k*n*a_start, 1);
    // C_1 = A_0 * A_0 // full karatsuba
    karatsuba(a_start, c_start+len3, 2*len2 + start, Layout, trans, n, k);
    // C_1 -= C_0
    cblas_daxpy(n*n, -1.0, c_double_array+n*n*c_start, 1, c_double_array+n*n*(c_start+len3), 1);
    // C_1 -= C_2
    cblas_daxpy(n*n, -1.0, c_double_array+n*n*(c_start+2*len3), 1, c_double_array+n*n*(c_start+len3), 1);
    // A_0 -= A_1 // clean up
    cblas_daxpy(k*n, -1.0, a_double_array+k*n*a_start, 1, a_double_array+k*n*(a_start+len2), 1);
  }
  
}

void mpmat::gradeschool(const int & a_start, const int & b_start, const int & c_start, const int & c_max,
                      CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE transa,
                      CBLAS_TRANSPOSE transb, const int & m, const int & n, const int & k,
                      const double alpha, const double beta){
  int start = a_start + b_start;
  int diff = c_max - start;
  //std::cout << "\tgradeschool start=(" << a_start << ", " << b_start << ") l=" << c_max << "\n";
  if(diff < 2){ // base case, just multiply
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
    //std::cout << "\t\tc = " << c_double_array[m*n*c_start] << "\n";
  }
  else { //if we don't need all four, then just do grade school
    int len2 = pow(2,ceil(log2(diff))-1);
    int len3 = pow(3,ceil(log2(diff))-1);
    // C_0 = A_0 * B_0 // karatsuba
    karatsuba(a_start, b_start, c_start, 2*len2 + start, Layout, transa, transb, m, n, k);
    // C_1 = A_1 * B_0 // grade school
    gradeschool(a_start+len2, b_start, c_start+len3, c_max, Layout, transa, transb, m, n, k);
    // C_2 = A_0 * B_1 // grade school // stored in C_2 temporarily!
    gradeschool(a_start, b_start+len2, c_start+2*len3, c_max, Layout, transa, transb, m, n, k);
    //std::cout << "\n";
    // C_1 += C_2
    cblas_daxpy(m*n*len3, 1.0, c_double_array+m*n*(c_start+2*len3), 1, c_double_array+m*n*(c_start+len3), 1);
    // C_2 = 0
    memset(c_double_array+m*n*(c_start+2*len3),0,m*n*len3*sizeof(double));
  }
}

void mpmat::gradeschool_bc(const int & a_start, const int & b_start, const int & c_start, const int & c_max,
                      CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE transa,
                      CBLAS_TRANSPOSE transb, const int & m, const int & n, const int & k,
                      const double alpha, const double beta){
  int start = a_start + b_start;
  int diff = c_max - start;
  //std::cout << "\tgradeschool_bc start=(" << a_start << ", " << b_start << ") l=" << c_max << "\n";
  if(diff < 2){ // base case, just multiply
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
    //std::cout << "\t\tc = " << c_double_array[m*n*c_start] << "\n";
  }
  else { //if we don't need all four, then just do grade school
    int len2 = pow(2,ceil(log2(diff))-1);
    int len3 = pow(3,ceil(log2(diff))-1);
    // C_0 = A_0 * B_0 // karatsuba
    karatsuba_bc(a_start, b_start, c_start, c_max, Layout, transa, transb, m, n, k);
    // C_1 = A_1 * B_0 // grade school
    gradeschool_bc(a_start+len2, b_start, c_start+len3, c_max, Layout, transa, transb, m, n, k);
    // C_2 = A_0 * B_1 // grade school // stored in C_2 temporarily!
    gradeschool_bc(a_start, b_start+len2, c_start+2*len3, c_max, Layout, transa, transb, m, n, k);
    //std::cout << "\n";
    // C_1 += C_2
    cblas_daxpy(m*n*len3, 1.0, c_double_array+m*n*(c_start+2*len3), 1, c_double_array+m*n*(c_start+len3), 1);
    // C_2 = 0
    memset(c_double_array+m*n*(c_start+2*len3),0,m*n*len3*sizeof(double));
  }
}

// IMPORTANT STEP: please make sure that a_double_array = b_double_array ahead of time
void mpmat::gradeschool(const int & a_start, const int & c_start, const int & c_max,
                      CBLAS_LAYOUT Layout, CBLAS_TRANSPOSE trans, const int & n, const int & k,
                      const double alpha, const double beta){ 
  int start = 2*a_start;
  int diff = c_max - start;
  if(diff < 2){ // base case, just multiply
   cblas_dsyrk(CblasColMajor,
                CblasUpper,(Layout == CblasRowMajor) != (trans == CblasTrans) ? CblasNoTrans : CblasTrans,
                n,k,alpha,
                a_double_array+k*n*a_start,
                (Layout == CblasRowMajor) != (trans == CblasTrans) ? n : k,
                beta,
                c_double_array+n*n*c_start,
                n);
  }
  else { //if we don't need all four, then just do grade school
    int len2 = pow(2,ceil(log2(diff))-1);
    int len3 = pow(3,ceil(log2(diff))-1);
    // C_0 = A_0 * B_0 // karatsuba
    karatsuba(a_start, c_start, c_max, Layout, trans, n, k);
    // C_1 = A_1 * B_0 // grade school
    gradeschool(a_start+len2, a_start, c_start+len3, c_max, Layout, trans, trans == CblasTrans ? CblasNoTrans : CblasTrans, n, n, k);
    // C_2 = A_0 * B_1 // grade school // stored in C_2 temporarily!
    gradeschool(a_start, a_start+len2, c_start+2*len3, c_max, Layout, trans, trans == CblasTrans ? CblasNoTrans : CblasTrans, n, n, k);
    // C_1 += C_2
    //cblas_daxpy(n*n, 1.0, c_double_array+n*n*(c_start+2*len3), 1, c_double_array+n*n*(c_start+len3), 1);
    // C_2 = 0
    memset(c_double_array+n*n*(c_start+2*len3),0,n*n*len3*sizeof(double));
  }
}

//given a length of a trinary tree, condense it into the nearly-binary tree that we need
//specifically a 3^n tree into a 2^(n+1)-1
void mpmat::treecondense(double * c, int size, int l){
  if (l < 9){
    return;
  }
  else {
    int len3 = l/9; // len3 = 3^(n-2), the length of the leaves two levels down.
                    // the len3 at the next level is 3*len3
    int n = log(l) / log(3) + .1; // tree level // the epsilon is because floating point sucks
    int len2 = (1 << (n-1)) -1; // the length of leaves two levels down post-squish
                            // the len2 at the next level is 2*len2 + 1

    //std::cout << "given l=" << l << ": len3=" << len3 << ", n=" << n << ", len2=" << len2 << "\n";

    // recurse! but carefully so as to not mess up the pointers
    treecondense(c+size*len3*6, size, l/3);
    treecondense(c+size*len3*3, size, l/3);
    treecondense(c, size, l/3);

    //if (l < 243){

     // now, a delicate operation: stitching together the newly squished answers.
    // Assuming that the previous level properly squished each answer,
    // the size went from 3^(n-1) to 2^n - 1. 
    // C_02 += C_10
    //std::cout << "\t\t\tadd array of length " << len2 << " from " << len3*3 << " to " << (len2+1) << "\n";
    cblas_daxpy(size*len2, 1.0, c+size*len3*3, 1, c+size*(len2+1), 1);
    // C_12 += C_20
    //std::cout << "\t\t\tadd array of length " << len2 << " from " << len3*6 << " to " << (len3*3+len2+1) << "\n";
    cblas_daxpy(size*len2, 1.0, c+size*len3*6, 1, c+size*(len3*3+len2+1), 1);

    // Fill in the gap between C_0 and C_1 without overlaps
    // i.e. overwrite the gap and where C_10 was
    //std::cout << "\t\t\tcopy array of length " << len2 << " from " << (len3*3+len2) << " to " << (2*len2+1) << "\n";
    cblas_dcopy(size*len2, c+size*(len3*3+len2), 1, c+size*(2*len2+1), 1);
    // Fill in the rest
    //std::cout << "\t\t\tcopy array of length " << 1 << " from " << (len3*3+2*len2) << " to " << (3*len2+1) << "\n";
    cblas_dcopy(size, c+size*(len3*3+2*len2), 1, c+size*(3*len2+1), 1);

    // get the last chunk. 3*len2+2 + len2+1 == 2^(n+1) - 1
    //std::cout << "\t\t\tcopy array of length " << len2+1 << " from " << len3*6+len2 << " to " << 3*len2+2 << "\n";
    cblas_dcopy(size*(len2+1), c+size*(len3*6+len2), 1, c+size*(3*len2+2), 1);
    
    // clean up a bit, juust in case
    //memset((c+size*((int)pow(2,n+1)-1)), 0, (l*size -pow(2,n+1)));
  //}
    return;
  }
  
}
