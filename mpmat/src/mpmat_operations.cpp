//
// Created by Petr Kravchuk on 8/14/17.
//

#include "mpmat.h"
#include <gmpxx.h>
#include <math.h>
#include <cassert>
#include <iostream>
#include "tests/mpmat_tests.h"
#include <mkl.h>
#include "Timers.h"
#include <cuda_runtime.h>
#include "cublas_v2.h"

template <typename T>
inline T ceil_div(T a, T b) {
    return a / b + ( (a % b != 0) & (a > 0) );
}

template <typename T>
inline T min(T a, T b) { return a < b ? a : b ; }

template <typename T>
inline T max (T a,T b) { return a>b ? a : b; }

void mpmatMultiplyGMPBaseCase(mpf_class & dest,
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
void cblas_dgemm_emulator(const int m,
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
}


void mpmat_gemm_reduced(
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

    std::cout << "Allocating double sizes " << mpmat_size_a << " " << mpmat_size_b << " " << mpmat_size_c << std::endl;
    std::cout << mpmat_size_a * m * k << std::endl;
    std::cout << mpmat_size_b * n * k << std::endl;
    std::cout << mpmat_size_c * m * n << std::endl;
    std::flush(std::cout);

    int mem_a = mpmat_size_a * m * k;
    int mem_b = mpmat_size_b * n * k;
    int mem_c = mpmat_size_c * m * n;

    auto a_double_array = new mpmat_double [mem_a];
    auto b_double_array = new mpmat_double [mem_b];

    auto c_double_array = new mpmat_double [mem_c];

    auto tmp            = new mpmat_double [ max( max(mem_a,mem_b), mem_c) ];

    memset(c_double_array, 0, mem_c * sizeof(mpmat_double));

    int expa, expb;

    std::cout << "Converting a to double" << std::endl;
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

    std::cout << "Converting b to double" << std::endl;
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

    std::cout << "Computing the product" << std::endl;
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

    std::cout << "Converting back" << std::endl;
    mpmatConvertDoubleToGMPVector(
            c,
            m*n,
            c_double_array,
            mpmat_size_c,
            mpmat_limb,
            expa+expb-mpmat_limb,
            tmp
    );

    delete [] a_double_array;
    delete [] b_double_array;
    delete [] c_double_array;
    delete [] tmp;

    timers["mpmat_gemm_reduced.complete"].stop();
}

void mpmat_gemm_reduced_gpu(
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

    std::cout << "Allocating double sizes " << mpmat_size_a << " " << mpmat_size_b << " " << mpmat_size_c << std::endl;
    std::cout << mpmat_size_a * m * k << std::endl;
    std::cout << mpmat_size_b * n * k << std::endl;
    std::cout << mpmat_size_c * m * n << std::endl;
    std::flush(std::cout);

    int mem_a = mpmat_size_a * m * k;
    int mem_b = mpmat_size_b * n * k;
    int mem_c = mpmat_size_c * m * n;

    auto a_double_array = new mpmat_double [mem_a];
    auto b_double_array = new mpmat_double [mem_b];

    auto c_double_array = new mpmat_double [mem_c];

    auto tmp            = new mpmat_double [ max( max(mem_a,mem_b), mem_c) ];

    memset(c_double_array, 0, mem_c * sizeof(mpmat_double));

    int expa, expb;

    std::cout << "Converting a to double" << std::endl;
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

    std::cout << "Converting b to double" << std::endl;
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

    
  int nDevices;

  cudaGetDeviceCount(&nDevices);
  for (int i = 0; i < nDevices; i++) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, i);
    printf("Device Number: %d\n", i);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n",
           prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n",
           prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
           2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
  }

    timers["mpmat_gemm_reduced.gpu_copy_forward"].start();

    mpmat_double * d_a, * d_b, *d_c, *h_a = new mpmat_double[mem_a], * h_b = new mpmat_double[mem_b];

    cudaMalloc(&d_a, mem_a*sizeof(mpmat_double));
    cudaMalloc(&d_b, mem_b*sizeof(mpmat_double));
    cudaMalloc(&d_c, mem_c*sizeof(mpmat_double));

    mpmat_double alpha = 1.0, beta = 1.0;

    cudaMemcpy(d_a,a_double_array,mem_a*sizeof(mpmat_double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_b,b_double_array,mem_b*sizeof(mpmat_double),cudaMemcpyHostToDevice);
    //cudaMemcpy(d_c,c_double_array,mem_c,cudaMemcpyHostToDevice);
    std::cout << "doing the GPU thing\n";

    cudaMemcpy(h_a,d_a,mem_a*sizeof(mpmat_double),cudaMemcpyDeviceToHost);
    cudaMemcpy(h_b,d_b,mem_b*sizeof(mpmat_double),cudaMemcpyDeviceToHost);

    for (int i = 0; i < mem_a; ++i)
      if(h_a[i] != a_double_array[i])
	std::cout << i << " in: " << a_double_array[i] << " out: " << h_a[i] << "\n";

    for (int i = 0; i < mem_b; ++i)
      if(h_b[i] != b_double_array[i])
	std::cout << i << " in: " << b_double_array[i] << " out: " << h_b[i] << "\n";
    

    timers["mpmat_gemm_reduced.gpu_copy_forward"].stop();

    cublasHandle_t handle;
    cublasCreate(&handle);

    timers["mpmat_gemm_reduced.multiplication"].start();

    //std::cout << "Computing the product" << std::endl;
    std::flush(std::cout);
    for (int i = 0; i < mpmat_size_c; i++) {
        for (int j = 0; j <= i; j++) {
	  cublasDgemm(handle,CUBLAS_OP_T,CUBLAS_OP_T,
		       m,
		       n,
		       k,
		       &alpha,
		       d_a+k*m*j,
		       Layout == CblasRowMajor ? k : m,
		       d_b+(i-j)*k*n,
		       Layout == CblasRowMajor ? n : k,
		       &beta,
		       d_c+i*m*n,
		       Layout == CblasRowMajor ? n : m
		       );
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

    timers["mpmat_gemm_reduced.gpu_copy_back"].start();

    mpmat_double * h_c = new mpmat_double[mem_c];

    // cudaMemcpy(h_a,d_a,mem_a,cudaMemcpyDeviceToHost);
    //cudaMemcpy(h_b,d_b,mem_b,cudaMemcpyDeviceToHost);
    cudaMemcpy(h_c,d_c,mem_c*sizeof(mpmat_double),cudaMemcpyDeviceToHost);

    /* cudaMemcpy(a_double_array,d_a,mem_a,cudaMemcpyDeviceToHost);
    cudaMemcpy(b_double_array,d_b,mem_b,cudaMemcpyDeviceToHost);
    cudaMemcpy(c_double_array,d_c,mem_c,cudaMemcpyDeviceToHost);*/

    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);

    timers["mpmat_gemm_reduced.gpu_copy_back"].stop();

    for (int i = 0; i < mem_c; ++i){
      if (h_c[i] != c_double_array[i])
	std::cout << "at position " << i << " expected " << c_double_array[i] << ", got " << h_c[i] << "\n";
    }

    delete [] h_c;

    cublasDestroy(handle);

    std::cout << "Converting back" << std::endl;
    mpmatConvertDoubleToGMPVector(
            c,
            m*n,
            c_double_array,
            mpmat_size_c,
            mpmat_limb,
            expa+expb-mpmat_limb,
            tmp
    );

    delete [] a_double_array;
    delete [] b_double_array;
    delete [] c_double_array;
    delete [] tmp;

    timers["mpmat_gemm_reduced.complete"].stop();
}

void mpmat_syrk_reduced(
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

  std::cout << "Allocating double sizes " << mpmat_size_a << " " << mpmat_size_c << std::endl;
  std::cout << mpmat_size_a * n * k << std::endl;
  std::cout << mpmat_size_c * n * n << std::endl;
  std::flush(std::cout);

  int mem_a = mpmat_size_a * n * k;
  int mem_c = mpmat_size_c * n * n;

  //std::cout << "mema is " << mema << " while memc is " << memc << "\n";

  auto a_double_array = new mpmat_double [mem_a];

  auto c_double_array = new mpmat_double [mem_c];

  auto tmp            = new mpmat_double [ max(mem_a, mem_c) ];

  int expa;
  
  std::cout << "Converting a to double" << std::endl;
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

  std::cout << "the exponent is " << expa << "\n";
  std::cout << "Computing the product" << std::endl;
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

  std::cout << "Converting back" << std::endl;
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
