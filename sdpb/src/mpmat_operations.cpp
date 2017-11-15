//
// Created by Petr Kravchuk on 8/14/17.
//
#define CUDA_API_PER_THREAD_DEFAULT_STREAM
#include "mpmat.h"
#include <gmpxx.h>
#include <math.h>
#include <cassert>
#include <iostream>
//#include "tests/mpmat_tests.h"
#include <mkl.h>
#include "Timers.h"
#include <omp.h>


#include <cuda_runtime.h>
#include "cublas_v2.h"
#include <cublasXt.h>

//mpmat_double *a_double_array,*b_double_array,*c_double_array,*tmp;

template <typename T>
inline T ceil_div(T a, T b) {
    return a / b + ( (a % b != 0) & (a > 0) );
}

template <typename T>
inline T min(T a, T b) { return a < b ? a : b ; }

template <typename T>
inline T max (T a,T b) { return a>b ? a : b; }

mpf_class * randomGMPVector(int size, int prec) {
    gmp_randclass rr(gmp_randinit_default);
    rr.seed(time(NULL));

    //mpf_class * tmp = static_cast<mpf_class*> ( malloc(size * sizeof(mpf_class)) );

    auto tmp = new mpf_class [size];

    for (int i = 0; i < size; i++) {
        //new(tmp+i) mpf_class("0",prec);
        tmp[i].set_prec(prec);
	//mpz_class z = rr.get_z_range(200);
	tmp[i] = 10*rr.get_f(prec)-5;
        //mpf_pow_ui(tmp[i].get_mpf_t(),tmp[i].get_mpf_t(), z.get_ui());
	
    }

    return tmp;
}


void mpmat::realloc(int mem_a, int mem_b, int mem_c){
  // std::cerr << "reallocing\n";
 if (mem_a > len_a){
   //std::cerr << "reallocing a\n";
      if (len_a != 0) {
	//	std::cerr << "len_a!=0\n";
	cudaFreeHost(a_double_array);
	//std::cerr << "deleted the relevant bits for a\n";
	for (int i = 0; i < gpu_count; ++i){
	  cudaSetDevice(i);
	  cudaFree(d_a[i]);
	  //std::cerr << "deleted d_a["<<i<<"]\n";
	}
      }
      //std::cerr << "allocing a\n";
      cudaMallocHost(&a_double_array,mem_a*sizeof(mpmat_double),cudaHostAllocPortable);
      //std::cerr << "alloc'd a\n";
      for (int i = 0; i < gpu_count; ++i){
	cudaSetDevice(i);
	cudaMalloc(d_a+i,mem_a*sizeof(mpmat_double));
      }
      len_a = mem_a;
    }
 //std::cerr << "realloc'd a\n";
 if (mem_b > len_b){

   //std::cerr << "reallocing b\n";
      if (len_b != 0) {
	cudaFreeHost(b_double_array);
	for (int i = 0; i < gpu_count; ++i){
	  cudaSetDevice(i);
	  cudaFree(d_b[i]);
	}
      }
      cudaMallocHost(&b_double_array,mem_b*sizeof(mpmat_double),cudaHostAllocPortable);
      for (int i = 0; i < gpu_count; ++i){
	cudaSetDevice(i);
	cudaMalloc(d_b+i,mem_b*sizeof(mpmat_double));
      }
      len_b = mem_b;
    }
 // std::cerr << "realloc'd b\n";
  if (mem_c > len_c){
    //std::cerr << "reallocing c\n";
      if (len_c != 0) {
	cudaFreeHost(c_double_array);
	for (int i = 0; i < gpu_count; ++i){
	  cudaSetDevice(i);
	  cudaFree(d_c[i]);
	}
      }
      cudaMallocHost(&c_double_array,mem_c*sizeof(mpmat_double),cudaHostAllocPortable);
      for (int i = 0; i < gpu_count; ++i){
	cudaSetDevice(i);
	cudaMalloc(d_c+i,mem_c*sizeof(mpmat_double));
      }
      len_c = mem_c;
    }
  //std::cerr << "realloc'd c\n";
    int mem_t = max(max(mem_a,mem_b),mem_c);
    if (mem_t > len_t){
      if (len_t != 0) delete [] tmp;
      tmp = new mpmat_double [mem_t];
      len_t = mem_t;
    }
}

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
    
#pragma omp parallel for
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

void mpmat::gemm_reduced_gpu(
        const CBLAS_LAYOUT Layout,
        const CBLAS_TRANSPOSE transa,
        const CBLAS_TRANSPOSE transb,
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

    timers["mpmat_gemm_reduced.complete"].resume();

    timers["mpmat_gemm_reduced.precalculations"].resume();

    int mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k)) )/ 2;
    int mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );
    int mpmat_size_b = ceil_div( abs(b[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );


    while ( 2 * mpmat_limb + ceil(log2(k*min(mpmat_size_a, mpmat_size_b))) > MPMAT_DOUBLE_MANT_IMPLICIT ) {
        mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k*min(mpmat_size_a, mpmat_size_b))) ) / 2;
        mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );
        mpmat_size_b = ceil_div( abs(b[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );
    }

    //std::cout << mpmat_limb << " is our limb for gemm\n";

    int mpmat_size_c = min(mpmat_size_a, mpmat_size_b);

    //std::cout << "Allocating double sizes " << mpmat_size_a << " " << mpmat_size_b << " " << mpmat_size_c << std::endl;
    //std::cout << mpmat_size_a * m * k << std::endl;
    //std::cout << mpmat_size_b * n * k << std::endl;
    //std::cout << mpmat_size_c * m * n << std::endl;
    std::flush(std::cout);

    int mem_a = mpmat_size_a * m * k;
    int mem_b = mpmat_size_b * n * k;
    int mem_c = mpmat_size_c * m * n;
    //int mem_t = max( max(mem_a,mem_b), mem_c);

    /*if (mem_a > len_a){
      if (len_a != 0) {
	delete [] a_double_array;
	cudaFree(d_a);
      }
      a_double_array = new mpmat_double [mem_a];
      cudaMalloc(&d_a,mem_a*sizeof(mpmat_double));
      len_a = mem_a;
    }
    if (mem_b > len_b){
      if (len_b != 0) {
	delete [] b_double_array;
	cudaFree(d_b);
      }
      b_double_array = new mpmat_double [mem_b];
      cudaMalloc(&d_b,mem_b*sizeof(mpmat_double));
      len_b = mem_b;
    }
    if (mem_c > len_c){
      if (len_c != 0) {
	delete [] c_double_array;
	cudaFree(d_c);
      }
      c_double_array = new mpmat_double [mem_c];
      cudaMalloc(&d_c,mem_c*sizeof(mpmat_double));
      len_c = mem_c;
    }
    if (mem_t > len_t){
      if (len_t != 0) delete [] tmp;
      tmp = new mpmat_double [mem_t];
      len_t = mem_t;
      }*/
    realloc(mem_a,mem_b,mem_c);

    memset(c_double_array, 0, mem_c * sizeof(mpmat_double));
    cudaMemset(d_c, 0, mem_c * sizeof(mpmat_double));

    timers["mpmat_gemm_reduced.precalculations"].stop();

    timers["mpmat_gemm_reduced.GMPtoDouble"].resume();

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

    timers["mpmat_gemm_reduced.GMPtoDouble"].stop();

    timers["mpmat_gemm_reduced.gpu_copy_forward"].resume();

    //mpmat_double * d_a, * d_b, *d_c;

    /*cudaMalloc(&d_a, mem_a*sizeof(mpmat_double));
    cudaMalloc(&d_b, mem_b*sizeof(mpmat_double));
    cudaMalloc(&d_c, mem_c*sizeof(mpmat_double));*/

    double alpha = 1.0, beta = 1.0;

    cudaMemcpy(d_a[0],a_double_array,mem_a*sizeof(mpmat_double),cudaMemcpyHostToDevice);
      cudaMemcpy(d_b[0],b_double_array,mem_b*sizeof(mpmat_double),cudaMemcpyHostToDevice);
    //cudaMemcpy(d_c,c_double_array,mem_c*sizeof(mpmat_double),cudaMemcpyHostToDevice);
    //std::cout << "doing the GPU thing\n";

    timers["mpmat_gemm_reduced.gpu_copy_forward"].stop();

    //cublasHandle_t handle;
    //cublasCreate(&handle);

    timers["mpmat_gemm_reduced.multiplication"].resume();

    //std::cout << "Computing the product" << std::endl;
    std::flush(std::cout);

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < mpmat_size_c; i++) {
        for (int j = 0; j <= i; j++) {
//cublasSetStream(handles[i],0);
	  //std::cout << "Thread: " << omp_get_thread_num() << "\n";
	  /* cublasXtDgemm(handle,
		      transa == CblasTrans ? CUBLAS_OP_N : CUBLAS_OP_T,
		      transb == CblasTrans ? CUBLAS_OP_N : CUBLAS_OP_T,
		       m,
		       n,
		       k,
		       &alpha,
		       a_double_array+k*m*j,
		       Layout == CblasRowMajor ? k : m,
		       b_double_array+(i-j)*k*n,
		       Layout == CblasRowMajor ? n : k,
		       &beta,
		       c_double_array+i*m*n,
		       Layout == CblasRowMajor ? n : m
		       );*/
	  cublasDgemm(handle,
		      (Layout == CblasRowMajor) != (transa == CblasTrans) ? CUBLAS_OP_T : CUBLAS_OP_N,
		      (Layout == CblasRowMajor) != (transb == CblasTrans) ? CUBLAS_OP_T : CUBLAS_OP_N,
		       m,
		       n,
		       k,
		       &alpha,
		       d_a[0]+k*m*j,
		      (Layout == CblasRowMajor) != (transa == CblasTrans) ? k : m,
						    d_b[0]+(i-j)*k*n,
		      (Layout == CblasRowMajor) != (transb == CblasTrans) ? n : k,
		       &beta,
		       d_c[0]+i*m*n,
		       Layout == CblasRowMajor ? n : m
		       );
	  /*cblas_dgemm(
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
		    );*/
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

    timers["mpmat_gemm_reduced.gpu_copy_back"].resume();

    //mpmat_double * h_c = new mpmat_double[mem_c];

    // cudaMemcpy(h_a,d_a,mem_a,cudaMemcpyDeviceToHost);
    //cudaMemcpy(h_b,d_b,mem_b,cudaMemcpyDeviceToHost);
    cudaMemcpy(c_double_array,d_c[0],mem_c*sizeof(mpmat_double),cudaMemcpyDeviceToHost);

    /* cudaMemcpy(a_double_array,d_a,mem_a,cudaMemcpyDeviceToHost);
    cudaMemcpy(b_double_array,d_b,mem_b,cudaMemcpyDeviceToHost);
    cudaMemcpy(c_double_array,d_c,mem_c,cudaMemcpyDeviceToHost);*/

    /*cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c)*/;

    timers["mpmat_gemm_reduced.gpu_copy_back"].stop();

    //for (int i = 0; i < mem_c; ++i){
    //if (h_c[i] != c_double_array[i])
    //	std::cout << "at position " << i << " expected " << c_double_array[i] << ", got " << h_c[i] << "\n";
    //}

    //delete [] h_c;

    //cublasDestroy(handle);

    timers["mpmat_gemm_reduced.DoubletoGMP"].resume();

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

    timers["mpmat_gemm_reduced.DoubletoGMP"].stop();

    /* delete [] a_double_array;
    delete [] b_double_array;
    delete [] c_double_array;
    delete [] tmp;*/

    timers["mpmat_gemm_reduced.complete"].stop();
}

/*void mpmat::syrk_reduced_gpu(
        const CBLAS_LAYOUT Layout,
        const CBLAS_TRANSPOSE transa,
        const int m,
        const int k,
        //const mpf_class alpha,
        const mpf_class * a,
        //const int lda,
        mpf_class * c
        //const int ldc
        ) {
  gemm_reduced_gpu(
        Layout,
        transa == CblasTrans ? CblasNoTrans : CblasTrans,
	transa,
        m,
	m,
        k,
        //const mpf_class alpha,
        a,
        //const int lda,
        a,
        //const int ldb,
        //const mpf_class beta,
        c
        //const int ldc
		   );
}*/

void mpmat::syrk_reduced_gpu(
        const CBLAS_LAYOUT Layout,
        const CBLAS_TRANSPOSE transa,
        const int m,
        const int k,
        //const mpf_class alpha,
        const mpf_class * a,
        //const int lda,
        mpf_class * c
        //const int ldc
        ) {

    timers["mpmat_syrk_reduced.complete"].resume();

    timers["mpmat_syrk_reduced.precalculations"].resume();

    int mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k)) )/ 2;
    int mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );


    while ( 2 * mpmat_limb + ceil(log2(k*mpmat_size_a)) > MPMAT_DOUBLE_MANT_IMPLICIT ) {
        mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k*mpmat_size_a)) ) / 2;
        mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );
    }

    int mpmat_size_c = mpmat_size_a;

    int mem_a = mpmat_size_a * m * k;
    int mem_c = mpmat_size_c * m * m;

    realloc(mem_a,max(mem_a,mem_c),mem_c);
    
    memset(c_double_array, 0, mem_c * sizeof(mpmat_double));
#pragma omp parallel for
    for (int i = 0; i < gpu_count; ++i){
      cudaSetDevice(i);
      cudaMemset(d_c[i], 0, mem_c * sizeof(mpmat_double));
    }

    timers["mpmat_syrk_reduced.precalculations"].stop();

    timers["mpmat_syrk_reduced.GMPtoDouble"].resume();

    int expa;
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

    timers["mpmat_syrk_reduced.GMPtoDouble"].stop();

    timers["mpmat_syrk_reduced.gpu_copy_forward"].resume();

    double alpha = 1.0, beta = 1.0;

#pragma omp parallel for 
    for (int i = 0; i < gpu_count; ++i){
      cudaSetDevice(i);
      cudaMemcpyAsync(d_a[i],a_double_array,mem_a*sizeof(mpmat_double),cudaMemcpyHostToDevice);
    }
    cudaThreadSynchronize();

    timers["mpmat_syrk_reduced.gpu_copy_forward"].stop();

    timers["mpmat_syrk_reduced.multiplication"].resume();
    std::flush(std::cout);
    

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < mpmat_size_c; i++) {
      int gpu_id = i * gpu_count / mpmat_size_c;
	cudaSetDevice(gpu_id);
      for (int j = 0; j < i / 2 + i % 2; j++) {
	
	  cublasDgemm(handles[gpu_id],
		      (Layout == CblasRowMajor) != (transa == CblasTrans) ? CUBLAS_OP_N : CUBLAS_OP_T,
		      (Layout == CblasRowMajor) != (transa == CblasTrans) ? CUBLAS_OP_T : CUBLAS_OP_N,
		       m,
		       m,
		       k,
		       &alpha,
		       d_a[gpu_id]+k*m*j,
		       (Layout == CblasRowMajor) != (transa == CblasTrans) ? m : k,
		       d_a[gpu_id]+(i-j)*k*m,
		       (Layout == CblasRowMajor) != (transa == CblasTrans) ? m : k,
		       &beta,
		       d_c[gpu_id]+i*m*m,
		       m
		       );
        }
	cublasDcopy(handles[gpu_id],m*m,d_c[gpu_id]+i*m*m,1,d_b[gpu_id]+i*m*m,1);
        cublasDgeam(handle,CUBLAS_OP_T,CUBLAS_OP_N,m,m,
		    &alpha,d_b[gpu_id]+i*m*m,m,
		    &beta,d_c[gpu_id]+i*m*m,m,
		    d_c[gpu_id]+i*m*m,m);
	// if significance of result is even, calculate the symmetric part
	if ( i % 2 == 0)
	  cublasDsyrk(handles[gpu_id],CUBLAS_FILL_MODE_UPPER,(Layout == CblasRowMajor) != (transa == CblasTrans) ? CUBLAS_OP_N : CUBLAS_OP_T,
		    m,k,
		    &alpha, d_a[gpu_id]+k*m*(i/2), (Layout == CblasRowMajor) != (transa == CblasTrans) ? m : k,
		    &beta, d_c[gpu_id]+m*m*i, m);
	 cudaMemcpyAsync(c_double_array+i*m*m,d_c[gpu_id]+i*m*m,m*m*sizeof(mpmat_double),cudaMemcpyDeviceToHost);
    }


    cudaThreadSynchronize();
    timers["mpmat_syrk_reduced.multiplication"].stop();

    timers["mpmat_syrk_reduced.gpu_copy_back"].resume();

    // gpu copy back is actually happening asynchronously with the multiplications
    timers["mpmat_syrk_reduced.gpu_copy_back"].stop();

    timers["mpmat_syrk_reduced.DoubletoGMP"].resume();

    mpmatConvertDoubleToGMPSymm(
            c,
            m,
            c_double_array,
            mpmat_size_c,
            mpmat_limb,
            expa+expa-mpmat_limb,
            tmp
				);

    timers["mpmat_syrk_reduced.DoubletoGMP"].stop();


    timers["mpmat_syrk_reduced.complete"].stop();
}

void mpmat::syrk_reduced(
        const CBLAS_LAYOUT Layout,
        const CBLAS_TRANSPOSE transa,
        const int m,
        const int k,
        //const mpf_class alpha,
        const mpf_class * a,
        //const int lda,
        mpf_class * c
        //const int ldc
        ) {

    //std::cout << "The non-c++ code is running, I swear\n";

    timers["mpmat_syrk_reduced.complete"].resume();

    timers["mpmat_syrk_reduced.precalculations"].resume();

    int mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k)) )/ 2;
    //int mpmat_limb = 1;
    int mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );


    while ( 2 * mpmat_limb + ceil(log2(k*mpmat_size_a)) > MPMAT_DOUBLE_MANT_IMPLICIT ) {
        mpmat_limb = ( MPMAT_DOUBLE_MANT_IMPLICIT - ceil(log2(k*mpmat_size_a)) ) / 2;
        mpmat_size_a = ceil_div( abs(a[0].get_mpf_t()->_mp_prec+1) * mp_bits_per_limb, mpmat_limb );
    }
    //std::cout << mpmat_limb << " is our limb for syrk\n";

    int mpmat_size_c = mpmat_size_a;

    //std::flush(std::cout);

    int mem_a = mpmat_size_a * m * k;
    int mem_c = mpmat_size_c * m * m;
    //std::cerr << "mem_a is " << mem_a << "\n";

    realloc(mem_a,max(mem_a,mem_c),mem_c);
    
    memset(c_double_array, 0, mem_c * sizeof(mpmat_double));
    //std::cout << "resetting d_c[i]\n";

    timers["mpmat_syrk_reduced.precalculations"].stop();

    timers["mpmat_syrk_reduced.GMPtoDouble"].resume();

    int expa;

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
  
    timers["mpmat_syrk_reduced.GMPtoDouble"].stop();

    timers["mpmat_syrk_reduced.gpu_copy_forward"].resume();

    double alpha = 1.0, beta = 1.0;

    //cudaMemcpy(d_c,c_double_array,mem_c*sizeof(mpmat_double),cudaMemcpyHostToDevice);
    //std::cout << "doing the GPU thing\n";

    timers["mpmat_syrk_reduced.gpu_copy_forward"].stop();

    //cublasHandle_t handle;
    //cublasCreate(&handle);

    timers["mpmat_syrk_reduced.multiplication"].resume();

    //std::cout << "Computing the product" << std::endl;
   
    

    for (int i = 0; i < mpmat_size_c; i++) {
      for (int j = 0; j < i / 2 + i % 2; j++) {
	cblas_dgemm(Layout,(transa == CblasTrans) ? CblasNoTrans : CblasTrans,transa,
		    m,m,k,alpha,a_double_array+k*m+j,(Layout == CblasRowMajor) != (transa == CblasTrans) ? m : k,
		    a_double_array+(i-j)*k*m,(Layout == CblasRowMajor) != (transa == CblasTrans) ? m : k,beta,
		    c_double_array+i*m*m,m);
	 
        }
      mkl_domatcopy('r','t',m,m,alpha,c_double_array+i*m*m,m,b_double_array+i*m*m,m);
      cblas_daxpy(m*m,alpha,b_double_array+i*m*m,1,c_double_array+i*m*m,1);
	//review this
	if ( i % 2 == 0)
	  cblas_dsyrk(Layout,CblasUpper,(Layout == CblasRowMajor) != (transa == CblasTrans) ? CblasTrans : CblasNoTrans, m, k, alpha, a_double_array+k*m*(i/2),(Layout == CblasRowMajor) != (transa == CblasTrans) ? m : k, beta, c_double_array+i*m*m,m);
    }


    //   std::cout << "done multiplying\n";
    timers["mpmat_syrk_reduced.multiplication"].stop();

    timers["mpmat_syrk_reduced.gpu_copy_back"].resume();

    //mpmat_double * h_c = new mpmat_double[mem_c];

    // cudaMemcpy(h_a,d_a,mem_a,cudaMemcpyDeviceToHost);
    //cudaMemcpy(h_b,d_b,mem_b,cudaMemcpyDeviceToHost);

    /*int gpu_imax[gpu_count+1];
    int gpu_id = -1;
    for (int i = 0; i < mpmat_size_c; ++i){
      int gpu = i * gpu_count / mpmat_size_c;
      if (gpu > gpu_id){
	gpu_id = gpu;
	gpu_imax[gpu] = i;
      }
    }
    gpu_imax[gpu_count] = mpmat_size_c;
#pragma omp parallel for shared(gpu_imax)
    for (int i = 0; i < gpu_count; ++i){
      cudaSetDevice(i);
      cudaMemcpy(c_double_array+gpu_imax[i]*m*m,d_c[i]+gpu_imax[i]*m*m,(gpu_imax[i+1]-gpu_imax[i])*m*m*sizeof(mpmat_double),cudaMemcpyDeviceToHost);
      }*/
    

    /* cudaMemcpy(a_double_array,d_a,mem_a,cudaMemcpyDeviceToHost);
    cudaMemcpy(b_double_array,d_b,mem_b,cudaMemcpyDeviceToHost);
    cudaMemcpy(c_double_array,d_c,mem_c,cudaMemcpyDeviceToHost);

    /*cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);*/

    timers["mpmat_syrk_reduced.gpu_copy_back"].stop();

    //for (int i = 0; i < mem_c; ++i){
    //if (h_c[i] != c_double_array[i])
    //	std::cout << "at position " << i << " expected " << c_double_array[i] << ", got " << h_c[i] << "\n";
    //}

    //delete [] h_c;

    //cublasDestroy(handle);

    timers["mpmat_syrk_reduced.DoubletoGMP"].resume();

    //std::cout << "Converting back" << std::endl;
    mpmatConvertDoubleToGMPSymm(
            c,
            m,
            c_double_array,
            mpmat_size_c,
            mpmat_limb,
            expa+expa-mpmat_limb,
            tmp
    );

    //for (int r = 0; r < m; ++r)
    //for (int col = 0; col < r; ++c)
    //	c[r*m + col] = c[col*m + r];

    timers["mpmat_syrk_reduced.DoubletoGMP"].stop();

    /* delete [] a_double_array;
    delete [] b_double_array;
    delete [] c_double_array;
    delete [] tmp;*/

    timers["mpmat_syrk_reduced.complete"].stop();
}
