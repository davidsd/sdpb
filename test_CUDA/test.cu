#include <sstream>
#include <iostream>
#include <stdio.h>
#include <cstdlib> 
#include <cuda_runtime_api.h>
#include <malloc.h>
#include <curand.h>
#include <cublas_v2.h>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include "omp.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp> 
#define uS_PER_SEC 1000000
#define uS_PER_mS 1000

typedef boost::numeric::ublas::matrix<double> matrix;

 
// Fill the array A(nr_rows_A, nr_cols_A) with random numbers on GPU

void GPU_fill_rand(double *A, int nr_rows_A, int nr_cols_A) {
     // Create a pseudo-random number generator
     curandGenerator_t prng;
     curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_DEFAULT);
 
     // Set the seed for the random number generator using the system clock
     curandSetPseudoRandomGeneratorSeed(prng, (unsigned long long) clock());
 
     // Fill the array with random numbers on the device
     curandGenerateUniformDouble(prng, A, nr_rows_A * nr_cols_A);
}




// Multiply the arrays A and B on GPU and save the result in C
// C(m,n) = A(m,k) * B(k,n)
// Also creates and destroys cuBLAS handle 
void gpu_blas_mmul(const double *A, const double *B, double *C, const int m, const int k, const int n) {
     int lda=m,ldb=k,ldc=m;
     const double alf = 1;
     const double bet = 0;
     const double *alpha = &alf;
     const double *beta = &bet;

     cublasHandle_t handle;
     cublasCreate(&handle);

      // Do the actual multiplication
     cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
     
     cublasDestroy(handle);
}



 // Multiply the arrays A and B on GPU and save the result in C
 // C(m,n) = A(m,k) * B(k,n)
void gpu_blas_mmul(const cublasHandle_t handle, const double *A, const double *B, double *C, const int m, const int k, const int n) {
     int lda=m,ldb=k,ldc=m;
     const double alf = 1;
     const double bet = 0;
     const double *alpha = &alf;
     const double *beta = &bet;
 
      // Do the actual multiplication
     cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}


// Multiply the array A with its transpse. 
// B(n, n) =  A(n, k) A^T(k, n)
void gpu_blas_mullWithTransp(const cublasHandle_t handle, const double *A, double *B, const int n, const int k) {
     int lda = n;
     const double alf = 1; 
     const double bet = 0;
     const double *alpha = &alf;
     const double *beta = &bet;

     
     // Do the actual multiplication
     cublasDsyrk(handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, n, k, alpha, A, lda, beta, B, lda);
}


//Print matrix A(nr_rows_A, nr_cols_A) storage in column-major format
void print_matrix(const double *A, int nr_rows_A, int nr_cols_A) {
 
     for(int i = 0; i < nr_rows_A; ++i){
         for(int j = 0; j < nr_cols_A; ++j){
             std::cout << A[j * nr_rows_A + i] << " ";
         }
         std::cout << std::endl;
     }
     std::cout << std::endl;
}

void moveMatrixToBoost(const double *A, int nr_rows_A, int nr_cols_A, matrix &B) {
  for(int i = 0; i < nr_rows_A; ++i){
    for(int j = 0; j < nr_cols_A; ++j){
      B (i, j) = A[j * nr_rows_A + i];
    }
  }
  std::cout << std::endl;
}



void print_memory() {
 // show memory usage of GPU

     size_t free_byte ;

     size_t total_byte ;

     cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;

     if ( cudaSuccess != cuda_status ){
            printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
            exit(1);
      }



     double free_db = (double)free_byte ;
     double total_db = (double)total_byte ;

     double used_db = total_db - free_db ;

     printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
                 used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

}


void matrixSquareIntoBlock(const double *A, double * B, const int nr_rows_A, const int nr_cols_A) {
  
  // #pragma omp parallel for schedule(dynamic)
#pragma omp parallel for schedule(dynamic)
  for (int c = 0; c < nr_rows_A; c++) {
    for (int r = 0; r <= c; r++) {
      double tmp = 0;
      for (int p = 0; p < nr_cols_A; p++)
        tmp += A[p * nr_rows_A + r] * A[p * nr_rows_A + c];
      B[r * nr_rows_A + c] = tmp;
      if (r != c)
        B[c * nr_rows_A + r] = tmp;
    }
  }
}

// Compute average difference for elements in matrix
double computeAverageDifference(const double *A, const double *B, const int nr_rows) {
  double sum = 0;
  for(int i = 0; i < nr_rows; ++i){
    for(int j = 0; j <= i; ++j){
      sum += abs(A[j * nr_rows + i] - B[j * nr_rows + i]);
    }
  }
  std::cout << "Average difference between the two is " << sum/(nr_rows * nr_rows) << std::endl;
  return(sum/(nr_rows * nr_rows));
}
 
int main(int argc, char *argv[]) {

    print_memory();
    int val1, val2; 

     if (argc >= 2)
     {
	std::istringstream iss1( argv[1] );
	std::istringstream iss2( argv[2] );
        if (iss1 >> val1)
        {
            // Conversion successful
        }
	if (iss2 >> val2)
	  {
            // Conversion successful                                                                                  
	  }
     }

     // Allocate 3 arrays on CPU     
     int nr_rows_A, nr_cols_A, nr_rows_B, nr_cols_B, nr_rows_C, nr_cols_C;
     timeval t1, t2;

     // for simplicity we are going to use square arrays
     nr_rows_A = nr_cols_B = nr_rows_C = nr_cols_C= val1;
     nr_cols_A = nr_rows_B = val2;
     
     double *h_A = (double *)malloc(nr_rows_A * nr_cols_A * sizeof(double));
     double *h_B = (double *)malloc(nr_rows_B * nr_cols_B * sizeof(double));
     double *h_C = (double *)malloc(nr_rows_C * nr_cols_C * sizeof(double));
     double *h_CNaive = (double *)malloc(nr_rows_C * nr_cols_C * sizeof(double));
     
     // Allocate 3 arrays on GPU
     double *d_A, *d_B, *d_C;
     cudaMalloc(&d_A,nr_rows_A * nr_cols_A * sizeof(double));
     cudaMalloc(&d_B,nr_rows_B * nr_cols_B * sizeof(double));
     cudaMalloc(&d_C,nr_rows_C * nr_cols_C * sizeof(double));
      
     // Fill the arrays A and B on GPU with random numbers
     std::cout << "Random filling..." << std::endl;
     GPU_fill_rand(d_A, nr_rows_A, nr_cols_A);
     GPU_fill_rand(d_B, nr_rows_B, nr_cols_B);
     std::cout << "Random filling done..." << std::endl;

     // Optionally we can copy the data back on CPU and print the arrays
     
     std::cout << "Transferring to memory..." << std::endl;
     // Starting the timer
     gettimeofday(&t1, NULL);

     cudaMemcpy(h_A,d_A,nr_rows_A * nr_cols_A * sizeof(double),cudaMemcpyDeviceToHost);
     cudaMemcpy(h_B,d_B,nr_rows_B * nr_cols_B * sizeof(double),cudaMemcpyDeviceToHost);
     //Ending the timer
     gettimeofday(&t2, NULL);
     float et1 = (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
     printf("Transfer GPU to memory time = %fms\n", et1);


     std::cout << "Transfer to memory ended..." << std::endl;

     std::cout << "A =" << std::endl;
     //print_matrix(h_A, nr_rows_A, nr_cols_A);
     std::cout << "B =" << std::endl;
     //print_matrix(h_B, nr_rows_B, nr_cols_B);
 


     // ***************************** Try matrix multiplication A A^T function that we've implemented 
     
     cublasHandle_t handle;
     cublasCreate(&handle);
     
     std::cout << "Multiplying matrix with transp..." << std::endl;
     // Starting the timer                                                                                           
     gettimeofday(&t1, NULL);

     gpu_blas_mullWithTransp(handle, d_A, d_C, nr_rows_A, nr_cols_A);

     //Ending the timer                                                                                              
     gettimeofday(&t2, NULL);
     et1 = (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
     printf("cuBLAS matrix multiplication time = %fms\n", et1);

     std::cout << "Multiplying matrix with transp done..." << std::endl;

     // ***************************** New multiplication requires handle to already be created 
     //std::cout << "Multiplying matrix..." << std::endl;
     //gpu_blas_mmul(handle, d_A, d_B, d_C, nr_rows_A, nr_cols_A, nr_cols_B);
     //std::cout << "Matrix multiplied..." << std::endl;
     cublasDestroy(handle);


     // ***************************** Old multiplication which also creates handle 
     // Multiply A and B on GPU
     //gpu_blas_mmul(d_A, d_B, d_C, nr_rows_A, nr_cols_A, nr_cols_B);
 
     // Copy (and print) the result on host memory
     std::cout << "Transferring matrix to memory...";
     cudaMemcpy(h_C,d_C,nr_rows_C * nr_cols_C * sizeof(double),cudaMemcpyDeviceToHost);
     std::cout << "Transferring matrix to memory ended...";
     std::cout << "C =" << std::endl;
     //print_matrix(h_C, nr_rows_C, nr_cols_C);
     std::cout << "Ended..." << std::endl;
     // ....

     print_memory();

     //Free GPU memory
     cudaFree(d_A);
     cudaFree(d_B);
     cudaFree(d_C);

     // ***************************** Check how long it would take to tranfer the matrices back into GPU memory
     
     cudaMalloc(&d_A,nr_rows_A * nr_cols_A * sizeof(double));
     // Starting the timer                                                                                           
     gettimeofday(&t1, NULL);
     cudaMemcpy(d_A,h_A,nr_rows_A * nr_cols_A * sizeof(float),cudaMemcpyHostToDevice);
     //Ending the timer                                                                                              
     gettimeofday(&t2, NULL);
     et1 = (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
     printf("Transfer CPU to GPU memory = %fms\n", et1);

     cudaFree(d_A);

     // ***************************** Naive multiplication that is currently implemented in SDPB                     
     // Starting the timer                                                                                           
     gettimeofday(&t1, NULL);

     std::cout << "Naive multiplication..." << std::endl;
     matrixSquareIntoBlock(h_A, h_CNaive, nr_rows_A, nr_cols_A);
     //Ending the timer                                                                                              
     gettimeofday(&t2, NULL);
     et1 = (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
     printf("Naive transpose multiplication time = %fms\n", et1);

     std::cout << "Naive multiplication ended..." << std::endl;
     std::cout << "C =" << std::endl;
     //print_matrix(h_C, nr_rows_C, nr_cols_C);            
     
     // Compare averages
     double diff = computeAverageDifference(h_C, h_CNaive, nr_rows_C);
     
     // **************************** Test boost as well 
     matrix m (nr_rows_A, nr_cols_A);
     matrix mT (nr_rows_A, nr_cols_A);
     matrix resultMatrix (nr_rows_A, nr_cols_A);
     std::cout << "Testing boost..." << std::endl;
     moveMatrixToBoost(h_A, nr_rows_A, nr_cols_A, m);
     mT = boost::numeric::ublas::trans(m);
     //std::cout << m << std::endl;
     std::cout << "Multiplication starting..." << std::endl;
     // Starting the timer                                                                                           
     gettimeofday(&t1, NULL);
     resultMatrix = boost::numeric::ublas::prod(m, mT);
     //Ending the timer                                                                                              
     gettimeofday(&t2, NULL);
     et1 = (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
     printf("Blas time = %fms\n", et1);

     std::cout << "C = " << std::endl;
     //std::cout << resultMatrix << std::endl;
     std::cout << "Multiplication ended..." <<std::endl;


     // Free memory
     free(h_A);
     free(h_B);
     free(h_C);
     free(h_CNaive);
     return 0;

 }
