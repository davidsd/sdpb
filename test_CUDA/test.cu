#include <sstream>
#include <iostream>
#include <stdio.h>
#include <cstdlib> 
#include <cuda_runtime_api.h>
#include <malloc.h>
#include <curand.h>
#include <cublas_v2.h>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include "omp.h"
#include <bitset>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp> 
#define uS_PER_SEC 1000000
#define uS_PER_mS 1000

typedef boost::numeric::ublas::matrix<double> matrix;


const int INT64L = 64;
const int DOUBLE_MANT = 53;

//Extracts up to 64 bits of one 64-bit variable (extractForm) starting from position1 up to position2                                                      
long long getBitsFromOneLong(const long long extractFrom, const int position1, const int position2) {
  assert(position1 <= position2);
  assert(INT64L >= position2);
  unsigned long long mask;
  if(position2 - position1 == INT64L) {mask = -1;} else {
    mask = ((1LL << (position2 - position1)) - 1) << (INT64L - position2);
  }
  return ((mask & extractFrom) >> (INT64L - position2));
}

//Extracts up to 64 bits of two 64-bits variables (extractFromA & extractFromB) starting from position1 up to position2                                    
long long getBitsFromTwoLong(const long long extractFromA, const long long extractFromB, const int position1, const int position2) {
  assert(position2 <= position1);
  return (getBitsFromOneLong(extractFromA, position1, INT64L) << position2)|(getBitsFromOneLong(extractFromB, 0, position2));
}

 
// Generates an array of 64-bit variables from mpf_class where each 64-bit variable have a maximum number of bits maxBits  
// The returned array of 64-bit variables is padded as if the exponent of the mpf_class variable is maxExp
void toBit(const long long a) {
  std::bitset<64> tmp(a);
  std::cout << tmp << std::endl;
}

// TODO: Eventually decouple allocation from calculation
//
void generateLongsFromGMP(const mpf_class a, long long *&x, int &sizeOfArray, const int ownLimbSize, const int padExp) {
  int size = a.get_mpf_t()->_mp_size; 
  // WARNING: _mp_exp is number of limbs, not actual exponent!
  // TODO: Test the fix below
  int realExp = a.get_mpf_t()->_mp_exp * INT64L;
  int padding = (padExp - realExp) / ownLimbSize;
  int padBitOffset = (padExp - realExp) % ownLimbSize;
  
  // Assert that the padding will work and that the maximumBitSize is given appropraitely 
  assert(realExp <= padExp);
  assert(ownLimbSize <= INT64L);

  // Find size of needed 64-bit variable array with padding and allocate it in memory
  sizeOfArray = (INT64L * size + (padExp - realExp)) / ownLimbSize;
  
  if ((INT64L * size + (padExp - realExp)) % ownLimbSize != 0) sizeOfArray += 1; 
  x = (long long *)malloc(sizeOfArray * sizeof(long long));
  
  // Add padding
  // Note that GMP stores the most significant digits of a number in _mp_d[size - 1]
  // That's why we start iterating the array of limbs in reverse
  for (int i = 0; i < padding; i++) x[i] = 0; 
  
  long long tmp  = a.get_mpf_t()->_mp_d[size - 1];
  x[padding] = getBitsFromOneLong(tmp, 0, ownLimbSize - padBitOffset);
  
  // Add all the elements in mpf_class to the result
  for (int i = padding + 1; i < sizeOfArray - 1; i++) {
    int leftGmpLimb  = size - 1 - (((i - padding) * ownLimbSize - padBitOffset) / INT64L);
    int leftBit      = ((i - padding) * ownLimbSize - padBitOffset) % INT64L;
    int rightGmpLimb = size - 1 - (((i - padding + 1) * ownLimbSize - padBitOffset) / INT64L);
    int rightBit     = ((i - padding + 1) * ownLimbSize - padBitOffset) % INT64L;
    // If true it means that all the bits are in the same limb. If flase it means that all the 
    // bits are in consecutive limbs.
    if (leftGmpLimb == rightGmpLimb) {
      long long tmp  = a.get_mpf_t()->_mp_d[leftGmpLimb];
      x[i] = getBitsFromOneLong(tmp, leftBit, rightBit);
    } else {
      long long tmpA = a.get_mpf_t()->_mp_d[leftGmpLimb];
      long long tmpB = a.get_mpf_t()->_mp_d[rightGmpLimb];
      x[i] = getBitsFromTwoLong(tmpA, tmpB, leftBit, rightBit);
    }
  }
  int leftBit = ((sizeOfArray - padding - 1) * ownLimbSize + padBitOffset) % INT64L;
  tmp = a.get_mpf_t()->_mp_d[0];
  x[sizeOfArray - 1] = getBitsFromOneLong(tmp, leftBit, INT64L);

  // TODO: Multiply longs by overall sign
}

// THIS HAS NOT BEEN TESTED YET
// TODO: decouple memory allocation from calculation eventually

void generateLongMatrixFromGMPMatrix(const mpf_class *a, const int nr_rows, const int nr_cols, long long **&x,  int &sizeOfArray, int &exp, const int ownLimbSize) {
  // Allocate memory for pointers that point to each element in matrix and to the array 
  // that gives the number of own limbs for each matrix element
  long long **tmpX;
  int * lengthOwnLimbs; 
  tmpX = (long long **) malloc( nr_rows * nr_cols * sizeof( long * ));
  lengthOwnLimbs = (int *) malloc( nr_rows * nr_cols * sizeof(int));

  // Find maximum exponent in matrix
  int maxExp = a[0].get_mpf_t()->_mp_exp;
  for(int i = 0; i < nr_rows; ++i)
    for(int j = 0; j < nr_cols; ++j) {
      int toCmp =  a[j * nr_rows + i].get_mpf_t()->_mp_exp;
      if (toCmp > maxExp) maxExp = toCmp;
    }
  exp = maxExp;
  
  // Generate the array of 64-bit matrices
  long minLengthOwnLimbs = LLONG_MAX;
  for(int i = 0; i < nr_rows; ++i)
    for(int j = 0; j < nr_cols; ++j) {
      generateLongsFromGMP(a[j * nr_rows + i], tmpX[j * nr_rows + i], lengthOwnLimbs[j * nr_rows + i], ownLimbSize, maxExp); 
      if (minLengthOwnLimbs > lengthOwnLimbs[j * nr_rows + i]) minLengthOwnLimbs = lengthOwnLimbs[j * nr_rows + i];
    }

  
  // Allocate the memory for the set of matrices (make all elements have the same length) such that 
  // elements of the matrices are closer in memory 
  // This might need to be rethough
  x = (long long **) malloc(minLengthOwnLimbs * sizeof *x + (minLengthOwnLimbs * (nr_rows * nr_cols * sizeof **x)));
  sizeOfArray = minLengthOwnLimbs;
  
  for(int k = 0; k < sizeOfArray; k++)
    for(int i = 0; i < nr_rows; ++i)
      for(int j = 0; j < nr_cols; ++j) {
	x[k][j * nr_rows + i] = tmpX[j * nr_rows + i][k];
      }
  free(tmpX);
  free(lengthOwnLimbs);
}
 
// This can be WAAAAAAY more optimized 
// Right now I think it takes time d^2 N^2. 
// It should take time N^2

// TODO: Add initialization functions for mpf_t's
// TODO: Maybe do bitwise operations and carry's by hand?
void addToGMP(mpf_class &a, const long long toAdd, const int bitToAdd) {
  mpf_t tmpDiv;
  mpf_init(tmpDiv);
  mpf_div_2exp(tmpDiv, a.get_mpf_t(), bitToAdd);
  mp_exp_t exp;
  std::cout << "Binary div..." << mpf_class(tmpDiv).get_str(exp, 2) << std::endl;
  std::cout << tmpDiv << std::endl;
  mpf_t tmpAdd;
  mpf_init(tmpAdd);
  mpf_add_ui(tmpAdd, tmpDiv, toAdd);
  std::cout << "Binary add..." << mpf_class(tmpAdd).get_str(exp, 2) << std::endl;
  mpf_t tmpMul;
  mpf_init(tmpMul);
  mpf_mul_2exp(tmpMul, tmpAdd, bitToAdd);
  std::cout << "Binary mul..." << mpf_class(tmpMul).get_str(exp, 2) << std::endl;
  a = mpf_class (tmpMul);
}

// NOT TESTED
void longToGMP(mpf_class &a, const long long *toAdd, const int size_toAdd, const int ownLimbSize, const int whereWeStart) {
  // Is the precision a global variable?... I assume so?
  a = mpf_class("0");
  for (int i = 0; i < size_toAdd; i++) {
    addToGMP(a, toAdd[i], whereWeStart - (i + 1) * ownLimbSize)
  }
}

// NOT TESTED
void numberMultiplicationBasecase(mpf_class &c, const mpf_class a, const mpf_class b) {
  if(a.get_prec() != b.get_prec()) {
    std::cout << "numberMultiplication::Numbers have different precision and therefore we will use the lower precision when multiplying" << std::endl;
  }
  
  long long *aS;
  int size_aS; 
  long long *bS; 
  int size_bS;
  
  // Define the maximum no of bits we can store in a 64-bit variable 
  // and the exponent to which we should pad one (or none) of the no
  int ownLimbSize = DOUBLE_MANT/2;
  int maxExp = max(a.get_mpf_t()->_m_exp, b.get_mpf_t()->_m_exp);
  generateLongsFromGMP(a, aS, size_aS, ownLimbSize, maxExp);
  generateLongsFromGMP(b, bS, size_bS, ownLimbSize, maxExp);
  int size = min(size_aS, size_bS);

  // Allocate memory to save the result in another array of 64-bit variables
  long long *res = (long long *)malloc(size * sizeof(long long)); 
  for (int i = 0; i < size; i++) {
    res[i] = 0;
    for (int j = 0; j < i + 1; j++) 
      res[i] += a[j] * b[i-j];
  }
}


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
     mpf_set_default_prec(100);
     mpf_class f("0.25");
     
     long long *mLimbs;
     int noOfmLimbs = 0;
     generateLongsFromGMP(f, mLimbs, noOfmLimbs, 7, 0);
     std::cout << "*****************************" << std::endl;
     mp_exp_t exp;
     std::cout << "Binary..." << f.get_str(exp, 2) << std::endl;
     std::cout << " Number of my limbs... " << noOfmLimbs << std::endl;
     for(int i = 0; i < noOfmLimbs; i++) {
       std::bitset<64> tr(mLimbs[i]);
       std::cout<< tr << std::endl;
     }
     std::cout << "*****************************" << std::endl;
     
     std::cout << "Binary..." << f.get_str(exp, 2) << std::endl;
     long long toAdd = 12;
     addToGMP(f, toAdd, -7);
     std::cout << "Binary added..." << f.get_str(exp, 2) << std::endl;
     std::bitset<64> tr(toAdd);
     std::cout << "what to add ..." << tr << std::endl; 

     std::cout << "*****************************" << std::endl;
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
     cudaThreadSynchronize();
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
     //matrixSquareIntoBlock(h_A, h_CNaive, nr_rows_A, nr_cols_A);
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
     //resultMatrix = boost::numeric::ublas::prod(m, mT);
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
