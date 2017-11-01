//
// Created by Petr Kravchuk on 8/14/17.
//

#include "mpmat.h"
#include <gmpxx.h>
#include <math.h>
#include <cassert>
#include <iostream>
//#include "tests/mpmat_tests.h"
#include <mkl.h>
#include "Timers.h"

//mpmat_double *a_double_array,*b_double_array,*c_double_array,*tmp;

template <typename T>
inline T ceil_div(T a, T b) {
    return a / b + ( (a % b != 0) & (a > 0) );
}

template <typename T>
inline T min(T a, T b) { return a < b ? a : b ; }

template <typename T>
inline T max (T a,T b) { return a>b ? a : b; }

void mpmat::mpmatMultiplyGMPBaseCase(mpf_class & dest,
                              const mpf_class  a,
                              const mpf_class  b) {

    int mpmat_limb = MPMAT_DOUBLE_MANT_IMPLICIT / 2;
    int mpmat_size_a = ceil_div( (a.get_mpf_t()->_mp_prec + 1) * mp_bits_per_limb, mpmat_limb );
    int mpmat_size_b = ceil_div( (b.get_mpf_t()->_mp_prec + 1) * mp_bits_per_limb, mpmat_limb );

    while ( 2 * mpmat_limb + ceil(log2(fmin(mpmat_size_a, mpmat_size_b))) > MPMAT_DOUBLE_MANT_IMPLICIT ) {
        mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(fmin(mpmat_size_a, mpmat_size_b))) ) / 2;
        mpmat_size_a = ceil_div( (a.get_mpf_t()->_mp_prec + 1) * mp_bits_per_limb, mpmat_limb );
        mpmat_size_b = ceil_div( (b.get_mpf_t()->_mp_prec + 1) * mp_bits_per_limb, mpmat_limb );
    }

    // Add some limbs for more reliable comparison with GMP. In principle we only know min(...) limbs of the result.
    int mpmat_size_c = min(mpmat_size_a, mpmat_size_b)+ceil_div(MPMAT_DOUBLE_MANT_IMPLICIT, mpmat_limb);

    assert(mpmat_limb > 0);

    mpmat_double * a_double = new mpmat_double[mpmat_size_a];
    mpmat_double * b_double = new mpmat_double[mpmat_size_b];
    mpmat_double * c_double = new mpmat_double[mpmat_size_c];

    int a_mp_exp = a.get_mpf_t() -> _mp_exp * mp_bits_per_limb;
    int b_mp_exp = b.get_mpf_t() -> _mp_exp * mp_bits_per_limb;

    mpmatConvertGMPToDouble(a,a_double,mpmat_size_a, mpmat_limb, a_mp_exp);
    mpmatConvertGMPToDouble(b,b_double,mpmat_size_b, mpmat_limb, b_mp_exp);

    for(int i = 0; i< mpmat_size_c; i++) {
        c_double[i]=0;
        for(int k = 0; k<=i; k++) {
            if (i-k < mpmat_size_a && k < mpmat_size_b) {
                c_double[i] += a_double[i - k] * b_double[k];
            }
        }
    }

    mpmatConvertDoubleToGMP(dest, c_double, mpmat_size_c,mpmat_limb, a_mp_exp+b_mp_exp-mpmat_limb);

    delete [] a_double;
    delete [] b_double;
    delete [] c_double;
}


// Test for row-major matrices, not optimzied
/*void cblas_dgemm_emulator(const int m,
                          const int n,
                          const int k,
                          const mpmat_double * a,
                          const mpmat_double * b,
                          mpmat_double * c){
    for (int ic = 0; ic < m; ic++) {
        for(int jc = 0; jc < n; jc++) {
            for(int l = 0; l < k; l++) {
                c[ n*ic + jc ] += a[ k*ic + l ] * b[ n*l + jc ];
            }
        }
    }
    }*/


void mpmat::gemm_reduced(
        const CBLAS_LAYOUT Layout,
        //const CBLAS_TRANSPOSE transa,
        //const CBLAS_TRANSPOSE transb,
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
        ) {
  ////std::cout << "The non-c++ code is running, I swear\n";

    timers["mpmat_gemm_reduced.complete"].start();

    int mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k)) )/ 2;
    int mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );
    int mpmat_size_b = ceil_div( abs(b[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );


    while ( 2 * mpmat_limb + ceil(log2(k*min(mpmat_size_a, mpmat_size_b))) > MPMAT_DOUBLE_MANT_IMPLICIT ) {
        mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k*min(mpmat_size_a, mpmat_size_b))) ) / 2;
        mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );
        mpmat_size_b = ceil_div( abs(b[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );
    }

    int mpmat_size_c = min(mpmat_size_a, mpmat_size_b);

    //std::cout << "Allocating double sizes " << mpmat_size_a << " " << mpmat_size_b << " " << mpmat_size_c << std::endl;
    //std::cout << mpmat_size_a * m * k << std::endl;
    //std::cout << mpmat_size_b * n * k << std::endl;
    //std::cout << mpmat_size_c * m * n << std::endl;
    std::flush(std::cout);

    int mem_a = mpmat_size_a * m * k;
    int mem_b = mpmat_size_b * n * k;
    int mem_c = mpmat_size_c * m * n;
    int mem_t = max( max(mem_a,mem_b), mem_c);

    if (mem_a > len_a){
      if (len_a != 0) delete [] a_double_array;
      a_double_array = new mpmat_double [mem_a];
      len_a = mem_a;
    }
    if (mem_b > len_b){
      if (len_b != 0) delete [] b_double_array;
      b_double_array = new mpmat_double [mem_b];
      len_b = mem_b;
    }
    if (mem_c > len_c){
      if (len_c != 0) delete [] c_double_array;
      c_double_array = new mpmat_double [mem_c];
      len_c = mem_c;
    }
    if (mem_t > len_t){
      if (len_t != 0) delete [] tmp;
      tmp = new mpmat_double [mem_t];
      len_t = mem_t;
    }

    memset(c_double_array, 0, mem_c * sizeof(mpmat_double));

    int expa, expb;

    //std::cout << "Converting a to double" << std::endl;
    std::flush(std::cout);
    mpmatConvertGMPToDoubleVector(
            a,
            m * k,
            a_double_array,
            mpmat_size_a,
            mpmat_limb,
            expa,
            tmp
    );

    //std::cout << "Converting b to double" << std::endl;
    std::flush(std::cout);
    mpmatConvertGMPToDoubleVector(
            b,
            n * k,
            b_double_array,
            mpmat_size_b,
            mpmat_limb,
            expb,
            tmp
    );

    timers["mpmat_gemm_reduced.multiplication"].start();

    //std::cout << "Computing the product" << std::endl;
    std::flush(std::cout);
    for (int i = 0; i < mpmat_size_c; i++) {
        for (int j = 0; j <= i; j++) {
            cblas_dgemm(
                    Layout,
                    CblasNoTrans,
                    CblasNoTrans,
                    m,
                    n,
                    k,
                    1,
                    a_double_array+k*m*j,
                    Layout == CblasRowMajor ? k : m,
                    b_double_array+(i-j)*k*n,
                    Layout == CblasRowMajor ? n : k,
                    1,
                    c_double_array+i*m*n,
                    Layout == CblasRowMajor ? n : m
            );
            /*cblas_dgemm_emulator(
                    m,
                    n,
                    k,
                    a_double_array+k*m*j,
                    b_double_array+(i-j)*k*n,
                    c_double_array+i*m*n
            );*/
        }
    }

    timers["mpmat_gemm_reduced.multiplication"].stop();

    //std::cout << "Converting back" << std::endl;
    mpmatConvertDoubleToGMPVector(
            c,
            m*n,
            c_double_array,
            mpmat_size_c,
            mpmat_limb,
            expa+expb-mpmat_limb,
            tmp
    );

    /* delete [] a_double_array;
    delete [] b_double_array;
    delete [] c_double_array;
    delete [] tmp;*/

    timers["mpmat_gemm_reduced.complete"].stop();
}

void mpmat::syrk_reduced(
			 const CBLAS_LAYOUT Layout,
			 const CBLAS_UPLO uplo,
			 const int n,
			 const int k,
			 const mpf_class * a,
			 mpf_class * c
			 ) {

  timers["mpmat_syrk_reduced.complete"].start();
  
  int mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k)) )/ 2;
  int mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );

  while ( 2 * mpmat_limb + ceil(log2(k*mpmat_size_a)) > MPMAT_DOUBLE_MANT_IMPLICIT ) {
    mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k*mpmat_size_a)) ) / 2;
    mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );
  }

  int mpmat_size_c = mpmat_size_a;

  //std::cout << "Allocating double sizes " << mpmat_size_a << " " << mpmat_size_c << std::endl;
  //std::cout << mpmat_size_a * n * k << std::endl;
  //std::cout << mpmat_size_c * n * n << std::endl;
  std::flush(std::cout);

  int mem_a = mpmat_size_a * n * k;
  int mem_c = mpmat_size_c * n * n;

  ////std::cout << "mema is " << mema << " while memc is " << memc << "\n";

  auto a_double_array = new mpmat_double [mem_a];

  auto c_double_array = new mpmat_double [mem_c];

  auto tmp            = new mpmat_double [ max(mem_a, mem_c) ];

  int expa;
  
  //std::cout << "Converting a to double" << std::endl;
  std::flush(std::cout);
  mpmatConvertGMPToDoubleVector(
			   a,
			   n * k,
			   a_double_array,
			   mpmat_size_a,
			   mpmat_limb,
			   expa,
			   tmp
			   );

  timers["mpmat_syrk_reduced.multiplication"].start();

  //std::cout << "the exponent is " << expa << "\n";
  //std::cout << "Computing the product" << std::endl;
  std::flush(std::cout);
  for (int i = 0; i < mpmat_size_c; i++) {
    for (int j = 0; j <= i; j++) {
      cblas_dsyrk(
		  Layout,
		  uplo,
		  CblasTrans,
		  n,
		  k,
		  1,
		  a_double_array+k*n*j,
		  Layout == CblasRowMajor ? k : n,
		  1,
		  c_double_array+i*n*n,
		  Layout == CblasRowMajor ? n : n
		  );
      /*cblas_dgemm_emulator(
	m,
	n,
	k,
	a_double_array+k*m*j,
	b_double_array+(i-j)*k*n,
	c_double_array+i*m*n
	);*/
    }
  }

  timers["mpmat_syrk_reduced.multiplication"].stop();

  //std::cout << "Converting back" << std::endl;
  mpmatConvertDoubleToGMPVector(
				c,
				n*n,
				c_double_array,
				mpmat_size_c,
				mpmat_limb,
				expa-mpmat_limb,
				tmp
				);
  timers["mpmat_syrk_reduced.complete"].stop();

  delete [] a_double_array;
  delete [] c_double_array;
  delete [] tmp;  
}
