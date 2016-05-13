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
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp> 
#define uS_PER_SEC 1000000
#define uS_PER_mS 1000

typedef boost::numeric::ublas::matrix<double> matrix;

using std::cout;
using std::endl;

const int INT64L = 64;
const int DOUBLE_MANT = 53;

int getSign(const mpf_class a) {
  if((a.get_mpf_t()->_mp_size) < 0) {return -1;}
  else {return 1;}
}

//Print matrix A(nr_rows_A, nr_cols_A) storage in column-major format                                                 
void print_matrix_long(const long long *A, int nr_rows_A, int nr_cols_A) {

  for(int i = 0; i < nr_rows_A; ++i){
    for(int j = 0; j < nr_cols_A; ++j){
      std::cout << A[j * nr_rows_A + i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void print_matrix_double(const double *A, int nr_rows_A, int nr_cols_A) {

  for(int i = 0; i < nr_rows_A; ++i){
    for(int j = 0; j < nr_cols_A; ++j){
      std::cout << (long long) A[j * nr_rows_A + i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

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
  int size = abs(a.get_mpf_t()->_mp_size); 
  // WARNING: _mp_exp is number of limbs, not actual exponent!
  // TODO: Test the fix below
  int realExp = a.get_mpf_t()->_mp_exp * INT64L;
  int padding = (padExp - realExp) / ownLimbSize;
  int padBitOffset = (padExp - realExp) % ownLimbSize;
  int sign  = getSign(a);

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
  x[padding] = sign * getBitsFromOneLong(tmp, 0, ownLimbSize - padBitOffset);
  
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
      x[i] = sign * getBitsFromOneLong(tmp, leftBit, rightBit);
    } else {
      long long tmpA = a.get_mpf_t()->_mp_d[leftGmpLimb];
      long long tmpB = a.get_mpf_t()->_mp_d[rightGmpLimb];
      x[i] = sign * getBitsFromTwoLong(tmpA, tmpB, leftBit, rightBit);
    }
  }
  int leftBit = ((sizeOfArray - padding - 1) * ownLimbSize + padBitOffset) % INT64L;
  tmp = a.get_mpf_t()->_mp_d[0];
  x[sizeOfArray - 1] = sign * getBitsFromOneLong(tmp, leftBit, INT64L);

  // TODO: Multiply longs by overall sign
}




// THIS HAS NOT BEEN TESTED YET
// TODO: decouple memory allocation from calculation eventually
void generateLongMatrixFromGMPMatrix(const mpf_class *a, long long **&x,int &sizeOfArray, const int nr_rows, const int nr_cols, int &exp, const int ownLimbSize) {
  // Allocate memory for pointers that point to each element in matrix and to the array 
  // that gives the number of own limbs for each matrix element
  long long **tmpX;
  int *lengthOwnLimbs; 
  tmpX = (long long **) malloc( nr_rows * nr_cols * sizeof( long * ));
  lengthOwnLimbs = (int *) malloc( nr_rows * nr_cols * sizeof(int));
  
  // Find maximum exponent in matrix
  int maxExp = a[0].get_mpf_t()->_mp_exp;
  for(int i = 0; i < nr_rows; ++i)
    for(int j = 0; j < nr_cols; ++j) {
      int toCmp =  a[j * nr_rows + i].get_mpf_t()->_mp_exp;
      if (toCmp > maxExp) maxExp = toCmp;
    }
  exp = maxExp * INT64L;
  
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
  sizeOfArray = minLengthOwnLimbs;
  x = (long long **)malloc(minLengthOwnLimbs * sizeof *x);
  for(int k = 0; k < sizeOfArray; k++)
    x[k] = (long long *) malloc(nr_rows * nr_cols * sizeof **x);
  
  for(int k = 0; k < sizeOfArray; k++)
    for(int i = 0; i < nr_rows; i++)
      for(int j = 0; j < nr_cols; j++) {
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
// Now tested: Seems to work!
mpf_class addToGMP(const mpf_class a, const long long toAdd, const int bitToAdd) {
  mpf_t tmpDiv;
  mpf_init(tmpDiv);
  if(bitToAdd < 0) {
    mpf_mul_2exp(tmpDiv, a.get_mpf_t(), abs(bitToAdd));
  } else {
    mpf_div_2exp(tmpDiv, a.get_mpf_t(), bitToAdd);
  }
  mpf_t tmpAdd;
  mpf_init(tmpAdd);
  mpf_add_ui(tmpAdd, tmpDiv, toAdd);
  mpf_t tmpMul;
  mpf_init(tmpMul);
  if (bitToAdd < 0) {
    mpf_div_2exp(tmpMul, tmpAdd, abs(bitToAdd));
  } else {
    mpf_mul_2exp(tmpMul, tmpAdd, bitToAdd);
  }
  return(mpf_class (tmpMul));
}


// Set a = a + b * 2^bitOffset. 
void addToMpf(mpf_t a, const long long b, const int bitOffset) {
  // bitOffset = limbOffset * GMP_NUMB_BITS + bitShift
  int limbOffset = bitOffset / GMP_NUMB_BITS;
  int bitShift   = bitOffset % GMP_NUMB_BITS;
  // ensure bitShift is positive
  if (bitShift < 0) {
    limbOffset -= 1;
    bitShift += GMP_NUMB_BITS;
  }

  unsigned long long bAbs = abs(b);

  // Let 2^GMP_NUMB_BITS = N. We would like to add/subtract
  // 
  //   bAbs * 2^bitOffset = (bAbs * 2^bitShift) * N^limbOffset
  //
  // So we write
  // 
  //   bAbs * 2^bitShift = head * 2^GMP_NUMB_BITS + tail
  unsigned long long head = bAbs >> (GMP_NUMB_BITS - bitShift);
  unsigned long long tail = bAbs << bitShift;

  // We now set
  //
  // a = ((a * N^(-limbOffset - 1) + head) * N + tail) * N^limbOffset
  //

  // a *= N^(-limbOffset - 1)
  a->_mp_exp -= limbOffset + 1;

  // a += head
  if (b > 0) {
    mpf_add_ui(a, a, head);
  } else {
    mpf_sub_ui(a, a, head);
  }

  // a *= N
  a->_mp_exp += 1;

  // a += tail
  if (b > 0) {
    mpf_add_ui(a, a, tail);
  } else {
    mpf_sub_ui(a, a, tail);
  }

  // a *= N^limbOffset
  a->_mp_exp += limbOffset;
}

void addToGMPMatrix(mpf_class *a, const long long *toAdd, const int nr_rows, const int nr_cols, const int bitToAdd) {
  #pragma omp parallel for schedule(dynamic)  
  for(int i = 0; i < nr_rows; ++i) {
    for(int j = 0; j < nr_cols; ++j) {
      addToMpf(a[j * nr_rows +  i].get_mpf_t(), toAdd[j * nr_rows +  i], bitToAdd);
      //a[j * nr_rows +  i] = addToGMP( a[j * nr_rows +  i], toAdd[j * nr_rows +  i], bitToAdd);
    }
  }
}

void addToGMPMatrixSymm(mpf_class *a, const long long *toAdd, const int nr_rows, const int bitToAdd) {
  #pragma omp parallel for schedule(dynamic)  
  for(int i = 0; i < nr_rows; ++i) {
    for(int j = 0; j <= i; ++j) {
      addToMpf(a[j * nr_rows +  i].get_mpf_t(), toAdd[j * nr_rows +  i], bitToAdd);
      //a[j * nr_rows +  i] = addToGMP( a[j * nr_rows +  i], toAdd[j * nr_rows +  i], bitToAdd);
    }
  }
}

// TO BE WRITTEN
void bitsToAddOneLong(unsigned long long &a, const unsigned long long b, 
		      const int bitToAdd, const int bitsToCopy) {
  unsigned long long mask;
  mask = ((1LL << (bitsToCopy)) - 1);
  a = ((mask & b) << bitToAdd);
}

// TO BE WRITTEN
void bitsToAddTwoLong(unsigned long long &a, unsigned long long &b, const unsigned long long c, 
		      const int bitToAdd, const int bitsToCopy) {
  unsigned long long mask1, mask2;
  //int ;
  //mask1 = ((1LL << ()) - 1);
  //mask2 = ((1LL << ()) - 1) << ();
  //a = ((mask2 & c) >> );
  //b = ((mask1 & c) << bitToAdd);
}

// TO BE WRITTEN
void addRemainder_InPlace(mpf_class &a, const long long remainder, const int limbGMPToAdd) {
  // if (limbGMPToAdd) {
  // }
  // if ( ) {
  //   a.get_mpf_t()->_mp_d[limbGMPToAdd] += remainder; 
  // } else {
  //   a.get_mpf_t()->_mp_d[limbGMPToAdd] += remainder; 
  //   addRemainder_InPlace(a, 1, limbGMPToAdd + 1);
  // }
}

// TO BE WRITTEN
void addToGMP_InPlace(mpf_class &a, const long long toAdd, const int bitToAdd) { 
  // Compare if the long that we have to add could actually exceed 
  int bitInLimb = bitToAdd % INT64L;
  int whichLimb = bitToAdd / (INT64L);
}

// NOT TESTED
void longToGMP(mpf_class &a, const long long *toAdd, const int size_toAdd, const int ownLimbSize, const int whereWeStart) {
  // Is the precision a global variable?... I assume so?
  a = mpf_class("0.0");
  for (int i = 0; i < size_toAdd; i++) {
    a = addToGMP(a, toAdd[i], whereWeStart - (i + 1) * ownLimbSize);
  }
}

// TESTED
void longToGMPMatrix(mpf_class *a, long long **toAdd, const int size_toAdd, const int nr_rows, const int nr_cols, const int ownLimbSize, const int whereWeStart) {
  for(int i = 0; i < nr_rows; ++i){
    for(int j = 0; j < nr_cols; ++j){
      a[j * nr_rows +  i] = mpf_class("0.0");
      for (int k = 0; k < size_toAdd; ++k) {
	a[j * nr_rows +  i]  = addToGMP(a[j * nr_rows +  i], toAdd[k][j * nr_rows +  i], whereWeStart - (k + 1) * ownLimbSize);
      }
    }
  }
} 

// TESTED                                                                                                                                                                                
void longToGMPMatrixDouble(mpf_class *a, double *toAdd, const int size_toAdd, const int nr_rows, const int nr_cols, const int ownLimbSize, const int whereWeStart) {
  for(int i = 0; i < nr_rows; ++i){
    for(int j = 0; j < nr_cols; ++j){
      a[j * nr_rows +  i] = mpf_class("0.0");
      for (int k = 0; k < size_toAdd; ++k) {
        a[j * nr_rows +  i]  = addToGMP(a[j * nr_rows +  i], (long long)toAdd[k * nr_rows * nr_cols + j * nr_rows +  i], whereWeStart - (k + 1) * ownLimbSize);
      }
    }
  }
}


// Now tested
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
  int maxExp = max(a.get_mpf_t()->_mp_exp, b.get_mpf_t()->_mp_exp) * INT64L;
  generateLongsFromGMP(a, aS, size_aS, ownLimbSize, maxExp);
  generateLongsFromGMP(b, bS, size_bS, ownLimbSize, maxExp);
  int size = min(size_aS, size_bS);
  
  // Allocate memory to save the result in another array of 64-bit variables
  long long *res = (long long *)malloc(size * sizeof(long long)); 
  for (int i = 0; i < size - 1; i++) {
    res[i] = 0;
    for (int j = 0; j < i + 1; j++) 
      res[i] += aS[j] * bS[i-j];
  }
  longToGMP(c, res, size - 1, ownLimbSize, 2 * maxExp - ownLimbSize);
  free(res);
}




void matrixProduct(long long *C, const long long *A, const long long *B, const int nr_rowsA, const int nr_colsA, const int nr_colsB) {

  // #pragma omp parallel for schedule(dynamic)                                                                                                                                 
  for (int c = 0; c < nr_rowsA; c++) {
    for (int r = 0; r < nr_colsB; r++) {
      long long tmp = 0;
      for (int p = 0; p < nr_colsA; p++)
        tmp += A[p * nr_rowsA + c] * B[r * nr_colsA + p];
      C[r * nr_rowsA + c] = tmp;
    }
  }
}

void matrixProductDouble(long long *C, const double *A, const double *B, const int nr_rowsA, const int nr_colsA, const int nr_colsB) {

  // #pragma omp parallel for schedule(dynamic)                                                                                                          
  for (int c = 0; c < nr_rowsA; c++) {
    for (int r = 0; r < nr_colsB; r++) {
      long long tmp = 0;
      for (int p = 0; p < nr_colsA; p++)
        tmp += (long long)(A[p * nr_rowsA + c] * B[r * nr_colsA + p]);
      C[r * nr_rowsA + c] = tmp;
    }
  }
}

void matrixProductGMP(mpf_class *C, const mpf_class *A, const mpf_class *B, const int nr_rowsA, const int nr_colsA, const int nr_colsB) {

  // #pragma omp parallel for schedule(dynamic)                                                                     
                                                                                                                    
  for (int c = 0; c < nr_rowsA; c++) {
    for (int r = 0; r < nr_colsB; r++) {
      mpf_class tmp = mpf_class("0");
      for (int p = 0; p < nr_colsA; p++)
        tmp += A[p * nr_rowsA + c] * B[r * nr_colsA + p];
      C[r * nr_rowsA + c] = tmp;
    }
  }
}


void matrixMultiplicationBasecase(mpf_class *c, mpf_class *a, mpf_class *b, long long **&aS, long long **&bS, const int nr_rowsA, const int nr_colsA, const int nr_colsB) {
  
  int size_aS = 0;
  int size_bS = 0;
  int exp = 0; 
  
  int ownLimbSize = DOUBLE_MANT/2 - ceil(log2((double) nr_colsA));
  
  generateLongMatrixFromGMPMatrix(a, aS, size_aS, nr_rowsA, nr_colsA, exp, ownLimbSize);  
  generateLongMatrixFromGMPMatrix(b, bS, size_bS, nr_colsA, nr_colsB, exp, ownLimbSize);
    
  int size = min(size_aS, size_bS);
  
  for(int i = 0; i < nr_rowsA; ++i){
    for(int j = 0; j < nr_colsB; ++j){
      c[j * nr_rowsA +  i] = mpf_class("0.0");
    }
  }
  
  
  long long *tmp = (long long *)malloc(nr_rowsA * nr_colsB * sizeof(long long));
  for (int i = 0; i < size - 1; i++) {
    for (int j = 0; j < i + 1; j++) {
      matrixProduct(tmp, aS[j], bS[i-j], nr_rowsA, nr_colsA, nr_colsB);
      addToGMPMatrix(c, tmp, nr_rowsA, nr_colsB, 2 * exp - (i + 2) * ownLimbSize); 
      }
  }
}


void generateRandomGMPMatrix(mpf_class *&a, const int nr_rows, const int nr_cols) {
  gmp_randclass rr(gmp_randinit_default);
  rr.seed(time(NULL));
  for(int i = 0; i < nr_rows; i++){
    for(int j = 0; j < nr_cols; j++){
      a[j * nr_rows +  i] = rr.get_f(500);
    }
  }
}

void printGMPMatrix(mpf_class *a, const int nr_rows, const int nr_cols) {
  std::cout << "Printing GMP matrix..." << std::endl;
  for(int i = 0; i < nr_rows; ++i){
    for(int j = 0; j < nr_cols; ++j){
      std::cout << a[j * nr_rows +  i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void printGMPMatrixDiff(mpf_class *a, mpf_class *b, const int nr_rows, const int nr_cols) {
  std::cout << "Printing GMP matrix difference..."<< std::endl;
  mp_exp_t exp;
  for(int i = 0; i < nr_rows; ++i){
    for(int j = 0; j < nr_cols; ++j){
      mpf_class tmp = (a[j * nr_rows +  i] - b[j * nr_rows +  i]);
      std::cout << tmp << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
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

// Fill the array A(nr_rows_A, nr_cols_A) with random numbers on GPU                                                  
void GPU_fill_rand_vec(unsigned long long *A, int length) {
  // Create a pseudo-random number generator                                                                       
  curandGenerator_t prng;
  curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_DEFAULT);

  // Set the seed for the random number generator using the system clock                                           
  curandSetPseudoRandomGeneratorSeed(prng, (unsigned long long) clock());

  // Fill the array with random numbers on the device                                                              
  curandGenerateLongLong(prng, A, length);
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
void gpu_blas_mulWithTransp(const cublasHandle_t handle, const double *A, double *B, const int n, const int k) {
  int lda = n;
  const double alf = 1; 
  const double bet = 0;
  const double *alpha = &alf;
  const double *beta = &bet;
  
     
  // Do the actual multiplication
  cublasDsyrk(handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, n, k, alpha, A, lda, beta, B, lda);
}


// Multiply the array A with its transpse. 
// B(n, n) =  A(n, k) A^T(k, n)
void gpu_blas_mulWithTranspAndSum(const cublasHandle_t handle, const double *A, const double *B, double *C, const int m, const int k) {
  int lda = n, ldb = k, ldc = k;
  const double alf = 1; 
  const double bet = 0;
  const double *alpha = &alf;
  const double *beta = &bet;

     
  // Do the actual multiplication
  cublasDsyrk(handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
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

__global__ void vecMult(long *a, long *b, long long res, int rem) {
  int id = blockIdx.x*blockDim.x+threadIdx.x;
  
}

// CUDA kernel. Each thread takes care of one element of c
__global__ void vecAdd(unsigned long long *a, unsigned long long *b, unsigned long long *c, int size)
{
  // Get our global thread ID
  int id = blockIdx.x*blockDim.x+threadIdx.x;
  
  // Make sure we do not go out of bounds
  if (id < size)
    c[id] = a[id] + b[id];
}


// CUDA kernel. Each thread takes care of one element of c                                                           
__global__ void vecAdd__wSign(double *a, long long *res, int size)
{
  // Get our global thread ID                                                                                       
  int id = blockIdx.x*blockDim.x+threadIdx.x;
  // Make sure we do not go out of bounds                                                                           
  if (id < size)
    res[id] += ((long long)a[id]);
}




// CUDA kernel. Each thread takes care of one element of c                                                          
__global__ void vecAdd_withRem(double *a, long long *res, long long *rem,  int size)
{
  // Get our global thread ID                                                                                       
  int id = blockIdx.x*blockDim.x+threadIdx.x;

  // Make sure we do not go out of bounds                                                                           
  if (id < size) {
    if (a[id] > 0 && res[id] > LLONG_MAX - a[id]) {
      /* handle overflow */
      rem[id] += 1;					       
      res[id] += (((long long)a[id]) - LLONG_MAX - 1);
    } else if (a[id] < 0 && res[id] < LLONG_MIN - a[id]) {
      /* handle underflow */
      rem[id] -= 1; 
      res[id] += (((long long)a[id]) - LLONG_MIN + 1);
    } else {
      /* handle regular cases */
      res[id] += a[id];
    }
  }
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


void estimateSize(const mpf_class *a, int &sizeOfArray, int &maxExp, const int nr_rows, const int nr_cols, const int ownLimbSize) {
  assert(ownLimbSize <= INT64L);

  int minSize = 0;
  maxExp = a[0].get_mpf_t()->_mp_exp;
  for(int i = 0; i < nr_rows; ++i)
    for(int j = 0; j < nr_cols; ++j) {
      int toCmp =  a[j * nr_rows + i].get_mpf_t()->_mp_exp;
      if (toCmp >= maxExp) {
        maxExp = toCmp;
        minSize = abs(a[j * nr_rows + i].get_mpf_t()->_mp_size);
      }
    }

  maxExp *= INT64L;
  sizeOfArray = (INT64L * minSize) / ownLimbSize;
  if ((INT64L * minSize) % ownLimbSize != 0) sizeOfArray += 1;
}



// All arrays need to already have allocated memory: 
// *a : array of mpf_class with size nr_rows * nr_cols
// *d_aS: array allocated in GPU with size sizeOfArray * nr_rows * nr_cols. Note that we flatten this 
//        array of matrices in order to speed up the access in the GPU.
//        Thus, to access the k-th limb from the (i, j) matrix entry one calls d_aS[k * nr_rows * nr_cols + j * nr_rows + i] 
// *tmpA: temporary array for doubles used for transfers that needs to be allocated nr_rows * nr_cols entries
// sizeOfArray: store number of own limbs that are saved for each entry
// nr_rows: number of matrix rows
// nr_cols: number of matrix columns 
// maxExp:  maximum power of 2 for the leading most limb among all the entries of the array *a 
void generateLongMatrixFromGMPMatrix_GPU(const mpf_class *a, double *d_aS, double *tmpA, 
					 int *sizeMatrix, int *realExpMatrix, int *signMatrix, int &sizeOfArray, 
					 const int nr_rows, const int nr_cols, int &maxExp, 
					 const int ownLimbSize) {
  assert(ownLimbSize <= INT64L);
  timeval t1, t2, tA, tB;
  double etAlloc = 0;
  double etAlgebra = 0;
  double etComputingSizes = 0;
  double etAddingLongs = 0;
  double etAddTrivial = 0;
  double etAccs = 0;
  #pragma omp parallel for schedule(dynamic)
  for(int j = 0; j < nr_cols; ++j) {
    for(int i = 0; i < nr_rows; ++i) {
      sizeMatrix[j * nr_rows + i] = abs(a[j * nr_rows + i].get_mpf_t()->_mp_size);
      realExpMatrix[j * nr_rows + i] = a[j * nr_rows + i].get_mpf_t()->_mp_exp * INT64L;
      signMatrix[j * nr_rows + i] =  getSign(a[j * nr_rows + i]);
    }
  }

  for (int k = 0; k < sizeOfArray - 1; k++) {
    gettimeofday(&t1, NULL);
    #pragma omp parallel for schedule(dynamic)
    for(int j = 0; j < nr_cols; ++j) {
      for(int i = 0; i < nr_rows; ++i) {
	//gettimeofday(&tA, NULL);
	int size = sizeMatrix[j * nr_rows + i];
	int realExp = realExpMatrix[j * nr_rows + i];
	int sign  = signMatrix[j * nr_rows + i];
	//gettimeofday(&tB, NULL);
	//etAccs += (((tB.tv_sec*uS_PER_SEC)+tB.tv_usec) - ((tA.tv_sec*uS_PER_SEC)+tA.tv_usec))/(float)uS_PER_mS;
	
	//gettimeofday(&tA, NULL);
	int padding = (maxExp - realExp) / ownLimbSize;
	int padBitOffset = (maxExp - realExp) % ownLimbSize;
	//gettimeofday(&tB, NULL);
	//etComputingSizes += (((tB.tv_sec*uS_PER_SEC)+tB.tv_usec) - ((tA.tv_sec*uS_PER_SEC)+tA.tv_usec))/(float)uS_PER_mSa;

	if(k > padding) {
	  //gettimeofday(&tA, NULL);
	  int leftGmpLimb  = size - 1 - (((k - padding) * ownLimbSize - padBitOffset) / INT64L);
	  int leftBit      = ((k - padding) * ownLimbSize - padBitOffset) % INT64L;
	  int rightGmpLimb = size - 1 - (((k - padding + 1) * ownLimbSize - padBitOffset) / INT64L);
	  int rightBit     = ((k - padding + 1) * ownLimbSize - padBitOffset) % INT64L;
	  //gettimeofday(&tB, NULL);
	  //etAddTrivial+= (((tB.tv_sec*uS_PER_SEC)+tB.tv_usec) - ((tA.tv_sec*uS_PER_SEC)+tA.tv_usec))/(float)uS_PER_mS;
	  //gettimeofday(&tA, NULL);
	  if (leftGmpLimb == rightGmpLimb) {
	    long long tmp  = a[j * nr_rows + i].get_mpf_t()->_mp_d[leftGmpLimb];
	    tmpA[j * nr_rows + i] = sign * getBitsFromOneLong(tmp, leftBit, rightBit);
	  } else {
	    long long tmp1 = a[j * nr_rows + i].get_mpf_t()->_mp_d[leftGmpLimb];
	    long long tmp2 = a[j * nr_rows + i].get_mpf_t()->_mp_d[rightGmpLimb];
	    tmpA[j * nr_rows + i] = sign * getBitsFromTwoLong(tmp1, tmp2, leftBit, rightBit);
	  }
	} else if (k < padding){
	  tmpA[j * nr_rows + i] = 0;
	} else {
	  long long tmp  = a[j * nr_rows + i].get_mpf_t()->_mp_d[size - 1];
	  tmpA[j * nr_rows + i] = sign * getBitsFromOneLong(tmp, 0, ownLimbSize - padBitOffset);
	}
	//gettimeofday(&tB, NULL);
	//etAddingLongs += (((tB.tv_sec*uS_PER_SEC)+tB.tv_usec) - ((tA.tv_sec*uS_PER_SEC)+tA.tv_usec))/(float)uS_PER_mS;
      } 
    }
    gettimeofday(&t2, NULL);
    etAlgebra +=  (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    // Transfer matrix to memeory
    //gettimeofday(&t1, NULL);
    cudaMemcpy(&d_aS[k * nr_rows * nr_cols], tmpA, nr_rows * nr_cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    //gettimeofday(&t2, NULL);
    //etAlloc +=  (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
  }

  #pragma omp parallel for schedule(dynamic)
  for(int i = 0; i < nr_rows; ++i) {
    for(int j = 0; j < nr_cols; ++j) {
      int realExp = a[j * nr_rows + i].get_mpf_t()->_mp_exp * INT64L;
      int padding = (maxExp - realExp) / ownLimbSize;
      int padBitOffset = (maxExp - realExp) % ownLimbSize;
      int sign  = getSign(a[j * nr_rows + i]);
      int leftBit = ((sizeOfArray - padding - 1) * ownLimbSize + padBitOffset) % INT64L;
      int tmp = a[j * nr_rows + i].get_mpf_t()->_mp_d[0];
      tmpA[j * nr_rows + i] = sign * getBitsFromOneLong(tmp, leftBit, INT64L);
    }
  }
  cudaMemcpy(&d_aS[(sizeOfArray - 1) * nr_cols * nr_rows], tmpA, nr_rows * nr_cols * sizeof(double),cudaMemcpyHostToDevice);
  //printf("Actual allocation time to GPU = %fms\n", etAlloc);
  printf("Algebra time in computing hands = %fms\n", etAlgebra);
  //printf("Calculating sizes = %fms\n", etComputingSizes);
  //printf("My function for obtaining the hands = %fms\n", etAddingLongs);
  //printf("Trivial stuff = %fms\n", etAddTrivial);
  //printf("Access time = %fms\n", etAccs);
}


// Estimates 
// maxFrac : 
// nr_rows :
// nr_cols :
int estimateMaxGPUAllocation(double maxFrac, int nr_rows, int nr_cols) {
  size_t free_byte ;
  size_t total_byte ;
  cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
  if ( cudaSuccess != cuda_status ){
    printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
    exit(1);
  }
  double free_db = ((double)free_byte) * maxFrac;
  return (int) (8 * free_db/(nr_rows * nr_cols * INT64L));
}


// Implements matrix multiplication c = a.b, and, according to whether or not the matrices fit 
// in GPU memory. If it does not fit we either apply the Karatsuba algorithm or the Toom-3 
// algorithm until all the matrices fit in GPU memory. 
void matrixMult_cuBlas(const cublasHandle_t handle, mpf_class *c, mpf_class *a, mpf_class *b, 
		       double *d_aS, double *d_bS, double *tmpTransferLongToGMP, 
		       long long *tmp, int *sizeMatrix, int *realExpMatrix, int *signMatrix,  
		       double *d_prodRes, long long *d_res, 
		       const int nr_rowsA, const int nr_colsA, const int nr_colsB) {
  int maxNoMatrices = estimateMaxGPUAllocation(0.8, (nr_rowsA + nr_colsB), nr_colsA);
  int size_aS = 0; 
  int size_bS = 0;
  int expA = 0;
  int expB = 0;
    
  int ownLimbSize = DOUBLE_MANT/2 - ceil(log2((double) nr_colsA) / 2);
  estimateSize(a, size_aS, expA, nr_rowsA, nr_colsA, ownLimbSize);
  estimateSize(b, size_bS, expB, nr_colsA, nr_colsB, ownLimbSize);
  
  if (max(size_aS, size_bS) + 1 < maxNoMatrices) {
    matrixMultiplicationBasecase_cuBlas(handle, c, a, b, d_aS, d_bS, tmpTransferLongToGMP, 
					tmp, sizeMatrix, realExpMatrix, signMatrix,  
					d_prodRes, d_res, nr_rowsA, nr_colsA, nr_colsB); 
  } else {
    toom2(handle, c, a, b, d_aS, d_bS, tmpTransferLongToGMP, 
	  tmp, sizeMatrix, realExpMatrix, signMatrix,  
	  d_prodRes, d_res, nr_rowsA, nr_colsA, nr_colsB);
  }
}

void matrixMultSymm_cuBlas(const cublasHandle_t handle, mpf_class *c, mpf_class *a, 
			   double *d_aS, double *tmpTransferLongToGMP, 
			   long long *tmp, int *sizeMatrix, int *realExpMatrix, int *signMatrix,  
			   double *d_prodRes, long long *d_res, 
			   const int nr_rowsA, const int nr_colsA, const float whatOrder) {
  int maxNoMatrices = estimateMaxGPUAllocation(0.8, nr_rowsA, nr_colsA);
  int size = 0; 
  int exp = 0;
 
  int ownLimbSize = DOUBLE_MANT/2 - ceil(log2((double) nr_colsA) / 2);
  
  estimateSize(a, size, exp, nr_rowsA, nr_colsA, ownLimbSize);

  if (size + 1 < maxNoMatrices) {
    matrixMultSymmBasecase_cuBlas(handle, c, a, d_aS, tmpTransferLongToGMP, 
				  tmp, sizeMatrix, realExpMatrix, signMatrix,  
				  d_prodRes, d_res, nr_rowsA, nr_colsA, prec); 
  } else {
    toom2(handle, c, a, b, d_aS, d_bS, tmpTransferLongToGMP, 
	  tmp, sizeMatrix, realExpMatrix, signMatrix,  
	  d_prodRes, d_res, nr_rowsA, nr_colsA, nr_colsB);
  }
}

// Implements Karatsuba algorithm for matrix multiplication c = a.b
// We split each matrix element in two parts which are equal in length up to the size of one GMP limb. 
void toom2(const cublasHandle_t handle, mpf_class *c, mpf_class *a, mpf_class *b, 
	   double *d_aS, double *d_bS, double *tmpTransferLongToGMP, 
	   long long *tmp, int *sizeMatrix, int *realExpMatrix, int *signMatrix,  
	   double *d_prodRes, long long *d_res, const int nr_rowsA, const int nr_colsA, const int nr_colsB) {
}

// Implements Toom-3 algorithm for matrix multiplciation c = a.b
// We split each matrix element in three parts which are equal in length up to the size of one GMP limb.
void toom3(const cublasHandle_t handle, mpf_class *c, mpf_class *a, mpf_class *b, 
	   double *d_aS, double *d_bS, double *tmpTransferLongToGMP, 
	   long long *tmp, int *sizeMatrix, int *realExpMatrix, int *signMatrix,  
	   double *d_prodRes, long long *d_res, const int nr_rowsA, const int nr_colsA, const int nr_colsB) {
}

// Implements Karatsuba algorithm for matrix multiplication c = a.a^T
//      
void toom2Symm(const cublasHandle_t handle, mpf_class *c, mpf_class *a, 
	       double *d_aS, double *tmpTransferLongToGMP, 
	       long long *tmp, int *sizeMatrix, int *realExpMatrix, int *signMatrix,  
	       double *d_prodRes, long long *d_res, const int nr_rowsA, const int nr_colsA, const float whatOrder) {
}

// Implements Toom-3 algorithm for matrix multiplication c = a.a^T
void toom3Symm(const cublasHandle_t handle, mpf_class *c, mpf_class *a, 
	       double *d_aS, double *tmpTransferLongToGMP, 
	       long long *tmp, int *sizeMatrix, int *realExpMatrix, int *signMatrix,  
	       double *d_prodRes, long long *d_res, const int nr_rowsA, const int nr_colsA, const float whatOrder) {
}

// Returns matrix product c = a.a^T where each entry is of the type mpf_class
//
// All arrays need to already have allocated memory:
// handle : Handle for cuBlas computations that needs to be previously allocated.                                                                   
// *c : array of mpf_class with size nr_rowsA * nr_rowsA                  
// *a : array of mpf_class with size nr_rowsA * nr_colsA
// *d_aS: array allocated in GPU with size sizeOfArray * nr_rowsA * nr_colsA. Note that we flatten this               
//        array of matrices in order to speed up the access in the GPU.                                              
//        Thus, to access the k-th limb from the (i, j) matrix entry one calls d_aS[k * nr_rowsA * nr_colsA + j * nr_rowsA + i]
// *tmpTransferLongToGMP : temporary array of doubles used for transfers that needs nr_rowsA * nr_colsA entries                    
// *tmp: temporary array of 64-bit vars used for transfers. Needs to be allocated nr_rowsA * nr_rowsA entries
// *d_prodRes: temporary matrix allocated in GPU memory which handles individual multiplications of matrices
// *d_res : temporary matrix allocated in GPU memory in which we sum up all individual multiplications of matrices
// *d_rem : in case we encounter overflow in the matrix d_res we use the entries in the matrix d_rem to store the remainders. 
//          Note that this matrix is only used when using the function vecAdd__wRem. When using vecAdd__wSign we assume that 
//          overflow never occurs. This is the case if the number of own generated limbs is greater than 1024. This correpons 
//          to a GMP precision ~20,000. For current bootstrap applications such precision is not needed.
// nr_rows_A: number of rows for matrix A as well as for matrix C                                                                                
// nr_cols_A: number of columns for matrix A as well as number of rows for matrix B
// nr_cols_B: number of columns for matrix B as well as for matrix C                                                                            
// maxExp:  maximum power of 2 for the leading most limb among all the entries of the array *a 
void matrixMultSymmBasecase_cuBlas(const cublasHandle_t handle, mpf_class *c, mpf_class *a, 
				   double *d_aS, double *tmpTransferLongToGMP, 
				   long long *tmp, int *sizeMatrix, int *realExpMatrix, int *signMatrix,  
				   double *d_prodRes, long long *d_res, const int nr_rowsA, const int nr_colsA, const float whatOrder) {
  
  int size = 0; 
  int expA = 0;
  timeval t1, t2;
  
  int ownLimbSize = DOUBLE_MANT/2 - ceil(log2((double) nr_colsA) / 2);
  gettimeofday(&t1, NULL);
  estimateSize(a, size_aS, expA, nr_rowsA, nr_colsA, ownLimbSize);
  
  generateLongMatrixFromGMPMatrix_GPU(a, d_aS, tmpTransferLongToGMP, sizeMatrix, realExpMatrix, signMatrix, 
				      size, nr_rowsA, nr_colsA, expA, ownLimbSize);
  gettimeofday(&t2, NULL);
  double etTransfer = (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
  
  std::cout << size << " " << expA << std::endl;

  for(int i = 0; i < nr_rowsA; ++i){
    for(int j = 0; j < nr_rowsA ++j){
      c[j * nr_rowsA +  i] = mpf_class("0.0");
      tmp[j * nr_rowsA +  i] = 0;
    }
  }

  // Number of threads in each thread block                                                                          
  int blockSize = 256; // Make sure this number is good? I don't know how
  // Number of thread blocks in grid                                                                                 
  int gridSize = (int)ceil((float)(nr_rowsA * nr_rowsA)/blockSize);
  
  double etMult = 0;
  double etAdd = 0;
  double etAlloc = 0;
  double etBackAlloc = 0;
  double etaddBack = 0;
  
  long sizeToCompute =  (whatOrder + maxExp) / ownLimbSize;
  
  for (int i = 0; i < min(size, sizeToCompute); i++) {
    gettimeofday(&t1, NULL);
    cudaMemcpy(d_res, tmp, nr_rowsA * nr_rowsA * sizeof(long long), cudaMemcpyHostToDevice);
    gettimeofday(&t2, NULL);
    etAlloc += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    //cudaMemcpy(d_rem, tmp, nr_rowsA * nr_colsB * sizeof(long long), cudaMemcpyHostToDevice);
    for (int j = 0; j < i + 1; j++) {
      gettimeofday(&t1, NULL);
      if (i - j != j)
	gpu_blas_mulWithTranspAndSum(handle, &d_aS[j * nr_rowsA * nr_rowsA], &d_aS[(i - j) * nr_colsA * nr_colsB], d_prodRes, nr_rowsA, nr_colsA); 
      else
	gpu_blas_mulWithTransp(handle, &d_aS[j * nr_rowsA * nr_rowsA], d_prodRes, nr_rowsA, nr_colsA); 
      cudaThreadSynchronize();
      gettimeofday(&t2, NULL);
      etMult += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
      
      gettimeofday(&t1, NULL);
      vecAdd__wSign<<<gridSize, blockSize>>>(d_prodRes, d_res,  nr_rowsA * nr_rowsA);
      cudaThreaMdSynchronize();
      gettimeofday(&t2, NULL);
      etAdd += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
      // This is safe with overflow until there are 1024 of our own limbs that need to be summed up
      // This case corresponds to a precision of ~20000 bits in GMP
    }
    gettimeofday(&t1, NULL);
    cudaMemcpy(tmp, d_res, nr_rowsA * nr_rowsA * sizeof(long long), cudaMemcpyDeviceToHost);
    gettimeofday(&t2, NULL);
    etBackAlloc += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    
    gettimeofday(&t1, NULL);
    addToGMPMatrixSymm(c, tmp, nr_rowsA, 2 * expA - (i + 2) * ownLimbSize);
    gettimeofday(&t2, NULL);
    etaddBack += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    
    #pragma omp parallel for schedule(dynamic)
    for(int l = 0; l < nr_rowsA; ++l){
      for(int k = 0; k < nr_rowsA; ++k){
	tmp[l *  nr_rowsA + k] = 0;
      }
    }
  }
  
   
  for (int i = size; i < min(2 * size - 1, sizeToCompute); i++) {
    gettimeofday(&t1, NULL);
    cudaMemcpy(d_res, tmp, nr_rowsA * nr_rowsA * sizeof(long long), cudaMemcpyHostToDevice);
    gettimeofday(&t2, NULL);
    etAlloc += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    //cudaMemcpy(d_rem, tmp, nr_rowsA * nr_colsB * sizeof(long long), cudaMemcpyHostToDevice);
    for (int j = i - size + 1; j < size; j++) {
      gettimeofday(&t1, NULL);
      if (i - j != j)
	gpu_blas_mulWithTranspAndSum(handle, &d_aS[j * nr_rowsA * nr_rowsA], &d_aS[(i - j) * nr_colsA * nr_colsB], d_prodRes, nr_rowsA, nr_colsA); 
      else
	gpu_blas_mulWithTransp(handle, &d_aS[j * nr_rowsA * nr_rowsA], d_prodRes, nr_rowsA, nr_colsA); 
      cudaThreadSynchronize();
      gettimeofday(&t2, NULL);
      etMult += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
      
      gettimeofday(&t1, NULL);
      vecAdd__wSign<<<gridSize, blockSize>>>(d_prodRes, d_res,  nr_rowsA * nr_rowsA);
      cudaThreaMdSynchronize();
      gettimeofday(&t2, NULL);
      etAdd += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
      // This is safe with overflow until there are 1024 of our own limbs that need to be summed up
      // This case corresponds to a precision of ~20000 bits in GMP
    }
    gettimeofday(&t1, NULL);
    cudaMemcpy(tmp, d_res, nr_rowsA * nr_rowsA * sizeof(long long), cudaMemcpyDeviceToHost);
    gettimeofday(&t2, NULL);
    etBackAlloc += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    
    gettimeofday(&t1, NULL);
    addToGMPMatrixSymm(c, tmp, nr_rowsA, 2 * expA - (i + 2) * ownLimbSize);
    gettimeofday(&t2, NULL);
    etaddBack += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    #pragma omp parallel for schedule(dynamic)
    for(int l = 0; l < nr_rowsA; ++l){
      for(int k = 0; k < nr_rowsA; ++k){
	tmp[l *  nr_rowsA + k] = 0;
      }
    }
  }
 
  printf("Transfer GPU = %fms\n", etTransfer);
  printf("Multiplication GPU = %fms\n", etMult);
  printf("Addition GPU = %fms\n", etAdd);
  printf("Alloc of zero to GPU = %fms\n", etAlloc);
  printf("Transfer from GPU to host = %fms\n", etBackAlloc);
  printf("Addition on CPU = %fms\n", etaddBack);
}




// Returns matrix product c = a.b where each entry is of the type mpf_class
//
// All arrays need to already have allocated memory:                                                                  
// *c : array of mpf_class with size nr_rowsA * nr_colsB                   
// *a : array of mpf_class with size nr_rowsA * nr_colsA
// *b : array of mpf_class with size nr_colsA * nr_colsB
// *d_aS: array allocated in GPU with size sizeOfArray * nr_rowsA * nr_colsA. Note that we flatten this               
//        array of matrices in order to speed up the access in the GPU.                                              
//        Thus, to access the k-th limb from the (i, j) matrix entry one calls d_aS[k * nr_rowsA * nr_colsA + j * nr_rowsA + i]
// *d_bS: array allocated in GPU with size sizeOfArray * nr_colsA * nr_colsB. Note that we flatten this                                         
//        array of matrices in order to speed up the access in the GPU.                                                                          
//        Thus, to access the k-th limb from the (i, j) matrix entry one calls d_bS[k * nr_colsA * nr_colsB + j * nr_colsA + i] 
// *tmpTransferLongToGMP : temporary array of doubles used for transfers that needs max(nr_rowsA * nr_colsA, nr_colsA * nr_colsB) entries        
// *tmp: temporary array of 64-bit vars used for transfers. Needs to be allocated nr_rowsA * nr_colsB entries
// *d_prodRes: temporary matrix allocated in GPU memory which handles individual multiplications of matrices
// *d_res : temporary matrix allocated in GPU memory in which we sum up all individual multiplications of matrices
// *d_rem : in case we encounter overflow in the matrix d_res we use the entries in the matrix d_rem to store the remainders. 
//          Note that this matrix is only used when using the function vecAdd__wRem. When using vecAdd__wSign we assume that 
//          overflow never occurs. This is the case if the number of own generated limbs is greater than 1024. This correpons 
//          to a GMP precision ~20,000. For current bootstrap applications such precision is not needed.
// nr_rows_A: number of rows for matrix A as well as for matrix C                                                                                
// nr_cols_A: number of columns for matrix A as well as number of rows for matrix B
// nr_cols_B: number of columns for matrix B as well as for matrix C                                                                             
// maxExp:  maximum power of 2 for the leading most limb among all the entries of the array *a 
void matrixMultiplicationBasecase_cuBlas(const cublasHandle_t handle, mpf_class *c, mpf_class *a, mpf_class *b, 
					 double *d_aS, double *d_bS, double *tmpTransferLongToGMP, 
					 long long *tmp, int *sizeMatrix, int *realExpMatrix, int *signMatrix,  
					 double *d_prodRes, long long *d_res,  //long long *d_rem, 
					 const int nr_rowsA, const int nr_colsA, const int nr_colsB) {
  
  int size_aS = 0; 
  int size_bS = 0;
  int expA = 0;
  int expB = 0;
  timeval t1, t2;
  
  int ownLimbSize = DOUBLE_MANT/2 - ceil(log2((double) nr_colsA) / 2);
  gettimeofday(&t1, NULL);
  estimateSize(a, size_aS, expA, nr_rowsA, nr_colsA, ownLimbSize);
  estimateSize(b, size_bS, expB, nr_colsA, nr_colsB, ownLimbSize);
  
  generateLongMatrixFromGMPMatrix_GPU(a, d_aS, tmpTransferLongToGMP, sizeMatrix, realExpMatrix, signMatrix, 
				      size_aS, nr_rowsA, nr_colsA, expA, ownLimbSize);
  generateLongMatrixFromGMPMatrix_GPU(b, d_bS, tmpTransferLongToGMP, sizeMatrix, realExpMatrix, signMatrix,
				      size_bS, nr_colsA, nr_colsB, expB, ownLimbSize);
  gettimeofday(&t2, NULL);
  double etTransfer = (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
  
  int size = min(size_aS, size_bS);
  std::cout << size_aS << " " << size_bS << " " << expA << " " << expB << std::endl;

  for(int i = 0; i < nr_rowsA; ++i){
    for(int j = 0; j < nr_colsB; ++j){
      c[j * nr_rowsA +  i] = mpf_class("0.0");
      tmp[j * nr_rowsA +  i] = 0;
    }
  }

  // Number of threads in each thread block                                                                          
  int blockSize = 256; // Make sure this number is good? I don't know how
  // Number of thread blocks in grid                                                                                 
  int gridSize = (int)ceil((float)(nr_rowsA * nr_colsB)/blockSize);
  
  double etMult = 0;
  double etAdd = 0;
  double etAlloc = 0;
  double etBackAlloc = 0;
  double etaddBack = 0;
  for (int i = 0; i < size; i++) {
    gettimeofday(&t1, NULL);
    cudaMemcpy(d_res, tmp, nr_rowsA * nr_colsB * sizeof(long long), cudaMemcpyHostToDevice);
    gettimeofday(&t2, NULL);
    etAlloc += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    //cudaMemcpy(d_rem, tmp, nr_rowsA * nr_colsB * sizeof(long long), cudaMemcpyHostToDevice);
    for (int j = 0; j < i + 1; j++) {
      gettimeofday(&t1, NULL);
      gpu_blas_mmul(handle, &d_aS[j * nr_rowsA * nr_colsA], &d_bS[(i - j) * nr_colsA * nr_colsB], d_prodRes, nr_rowsA, nr_colsA, nr_colsB);
      cudaThreadSynchronize();
      gettimeofday(&t2, NULL);
      etMult += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
      
      gettimeofday(&t1, NULL);
      vecAdd__wSign<<<gridSize, blockSize>>>(d_prodRes, d_res,  nr_rowsA * nr_colsB);
      cudaThreaMdSynchronize();
      gettimeofday(&t2, NULL);
      etAdd += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
      // This is safe with overflow until there are 1024 of our own limbs that need to be summed up
      // This case corresponds to a precision of ~20000 bits in GMP
    }
    gettimeofday(&t1, NULL);
    cudaMemcpy(tmp, d_res, nr_rowsA * nr_colsB * sizeof(long long), cudaMemcpyDeviceToHost);
    gettimeofday(&t2, NULL);
    etBackAlloc += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    
    gettimeofday(&t1, NULL);
    addToGMPMatrix(c, tmp, nr_rowsA, nr_colsB, expA + expB - (i + 2) * ownLimbSize);
    gettimeofday(&t2, NULL);
    etaddBack += (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    //cudaMemcpy(tmp, d_rem, nr_rowsA * nr_colsB * sizeof(long long),cudaMemcpyDeviceToHost);
    //addToGMPMatrix(c, tmp, nr_rowsA, nr_colsB, 2 * exp - (i + 2) * ownLimbSize + (INT64L - 1));
    #pragma omp parallel for schedule(dynamic)
    for(int k = 0; k < nr_rowsA; ++k){
      for(int l = 0; l < nr_colsB; ++l){
	tmp[l *  nr_rowsA + k] = 0;
      }
    }
  }
  printf("Transfer GPU = %fms\n", etTransfer);
  printf("Multiplication GPU = %fms\n", etMult);
  printf("Addition GPU = %fms\n", etAdd);
  printf("Alloc of zero to GPU = %fms\n", etAlloc);
  printf("Transfer from GPU to host = %fms\n", etBackAlloc);
  printf("Addition on CPU = %fms\n", etaddBack);
}


int estimateMaxGPUAllocation(double maxFrac, int nr_rows, int nr_cols) {
  size_t free_byte ;
  size_t total_byte ;
  cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
  if ( cudaSuccess != cuda_status ){
    printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
    exit(1);
  }
  double free_db = ((double)free_byte) * maxFrac;
  return (int) (8 * free_db/(nr_rows * nr_cols * INT64L));
}


int testAddToMpf() {
  mpf_set_default_prec(300);
  mpf_class x("3.14159265358");
  long long a = 12345;
  long long b = -4321899;
  cout.precision(300);
  cout << x << endl;
  addToMpf(x.get_mpf_t(), a, 62);
  cout << x << endl;
  addToMpf(x.get_mpf_t(), b, -4);
  cout << x << endl;
}
 
int lucaGiantTest(int argc, char *argv[]) {
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
     mpf_set_default_prec(300);
     
     timeval t1, t2;

     mpf_class f("-3.23124");
     mpf_class fCopy = f;
     mp_exp_t exp;
     long long *mLimbs;
     int noOfmLimbs = 0;
     //std::cout << getSign(f) << std::endl;
     //std::cout << "Res..." << f.get_str(exp, 10) << std::endl;
     //generateLongsFromGMP(f, mLimbs, noOfmLimbs, 7, 10);
     std::cout << "*****************************" << std::endl;
     //std::cout << "Binary..." << f.get_str(exp, 2) << std::endl;
     //std::cout << " Number of my limbs... " << noOfmLimbs << std::endl;
     //for(int i = 0; i < noOfmLimbs; i++) {
     //  std::bitset<64> tr(mLimbs[i]);
     //  std::cout<< tr << std::endl;
     //}
     std::cout << "*****************************" << std::endl;
     mpf_class getResBack;
     //longToGMP(getResBack, mLimbs, noOfmLimbs, 7, 10);
     //std::cout << getResBack.get_str(exp, 2) << std::endl;
     //std::cout << getResBack - f << std::endl;
     //std::cout << "*****************************" << std::endl;
     
     //std::cout << "Binary..." << f.get_str(exp, 2) << std::endl;
     //long long toAdd = 11232145214524552;
     //mpf_class result = addToGMP(f, toAdd, -5);
     //std::cout << "Binary added..." << result.get_str(exp, 2) << std::endl;
     //std::bitset<64> tr(toAdd);
     //std::cout << "what to add ..." << tr << std::endl; 
     //std::cout << "Binary added..." << result.get_str(exp, 10) << std::endl;
     std::cout << "*****************************" << std::endl;
     //std::cout << "Testing...";
     //mpf_class a1("0.23124");
     //mpf_class a2("0.251253124");
     //mpf_class a3; 

     //numberMultiplicationBasecase(a3, a2, a1);
     //std::cout << a3 << std::endl;
     //std::cout << a3 - (a1 * a2) << std::endl;
     std::cout << "*****************************" << std::endl;
     int nr_rows_A, nr_cols_A, nr_rows_B, nr_cols_B, nr_rows_C, nr_cols_C;
     std::cout << "Generating random matrix" << std::endl;
     nr_rows_A = nr_cols_B = nr_rows_C = nr_cols_C= val1;
     nr_cols_A = nr_rows_B = val2;
     std::cout << "Allocating..." << nr_rows_A << " " << nr_cols_A <<  std::endl;
     std::cout << sizeof(mpf_class) << std::endl;
     mpf_class *randA = new mpf_class[nr_rows_A * nr_cols_A];
     mpf_class *randB = new mpf_class[nr_rows_B * nr_cols_B];
     mpf_class *randC = new mpf_class[nr_rows_C * nr_cols_C];
     mpf_class *randCNaive = new mpf_class[nr_rows_C * nr_cols_C];
     mpf_class *randCGPU = new mpf_class[nr_rows_C * nr_cols_C];
     std::cout << "Allocated..." << std::endl;
     generateRandomGMPMatrix(randA, nr_rows_A, nr_cols_A);
     generateRandomGMPMatrix(randB, nr_rows_B, nr_cols_B);
     generateRandomGMPMatrix(randC, nr_rows_C, nr_cols_C);
     //printGMPMatrix(randA, nr_rows_A, nr_cols_A);
     std::cout << "*****************************" << std::endl;
     long long ** GMPtoLong;
     long long ** GMPtoLong2;
     int maxExpo = 0;
     //std::cout << "Generating an array of longs" << std::endl;
     //generateLongMatrixFromGMPMatrix(randA, GMPtoLong, noOfmLimbs, nr_rows_A, nr_cols_A, maxExpo, 7);
     //generateLongMatrixFromGMPMatrix(randA, GMPtoLong2, noOfmLimbs, nr_rows_A, nr_cols_A, maxExpo, 7);
     //std::cout << "Generated..." << std::endl;
     //std::cout << noOfmLimbs << std::endl;
     //std::cout << maxExpo << std::endl;
     std::cout << "*****************************" << std::endl;
     mpf_class *randACopy = new mpf_class[nr_rows_A * nr_cols_A];
     longToGMPMatrix(randACopy, GMPtoLong, noOfmLimbs, nr_rows_A, nr_cols_A, 7, maxExpo);
     //printGMPMatrix(randACopy, nr_rows_A, nr_cols_A);
     //std::cout << randACopy[0].get_prec() << std::endl;
     //std::cout << randACopy[0].get_str(exp, 10) << std::endl;
     //std::cout << randA[0].get_str(exp, 10) << std::endl;
     //printGMPMatrixDiff(randACopy, randA, nr_rows_A, nr_cols_A);
     std::cout << "*****************************" << std::endl;
     std::cout << "Testing matrix multiplication No GPU..." << std::endl;
     long long ** GMPtoLongA;
     long long ** GMPtoLongB;
     //printGMPMatrix(randA, nr_rows_A, nr_cols_A);
     //printGMPMatrix(randB, nr_cols_A, nr_cols_B);
     //matrixMultiplicationBasecase(randC, randA, randB, GMPtoLongA, GMPtoLongB, nr_rows_A, nr_cols_A, nr_cols_B);
     //printGMPMatrix(randC, nr_rows_C, nr_cols_C);
     gettimeofday(&t1, NULL);
     matrixProductGMP(randCNaive, randA, randB, nr_rows_A, nr_cols_A, nr_cols_B);
     gettimeofday(&t2, NULL);
     double etCPU = (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
     //printGMPMatrixDiff(randC, randCNaive, nr_rows_C, nr_cols_C);
     std::cout << "*****************************" << std::endl;
     std::cout << "Testing vector addition on GPU" << std::endl;
     long long *a  = (long long *)malloc(nr_rows_A * nr_cols_A * sizeof(long long));
     long long *b = (long long *)malloc(nr_rows_A * nr_cols_A * sizeof(long long));
     unsigned long long *d_a; 
     unsigned long long *d_b; 
     unsigned long long *d_c;
     //cudaMalloc(&d_a, nr_rows_A * nr_cols_A * sizeof(long long));
     //cudaMalloc(&d_b, nr_rows_A * nr_cols_A * sizeof(long long));
     //cudaMalloc(&d_c, nr_rows_A * nr_cols_A * sizeof(long long));
     //GPU_fill_rand_vec(d_a, nr_rows_A * nr_cols_A);
     //GPU_fill_rand_vec(d_b, nr_rows_A * nr_cols_A);
     // Number of threads in each thread block
     //int blockSize = 1024;
     // Number of thread blocks in grid
     //int gridSize = (int)ceil((float)(nr_rows_A * nr_cols_A)/blockSize);
     
     //vecAdd<<<gridSize, blockSize>>>(d_a, d_b, d_c, nr_rows_A * nr_cols_A);
     //vecAdd_withRem<<<gridSize, blockSize>>>(d_c, d_a, nr_rows_A * nr_cols_A);
     long long tryNegative = 1;
     double tryNeg = -1;
     double tryNeg2 = -2; 
     //std::bitset<64> tryN(tryNeg * tryNeg2);
     //std::cout << tryN << std::endl;
     
     std::cout << "*****************************" << std::endl;
     std::cout << "*** TESTING GMP MATRIX TO LONG MATRIX ON GPU ***" << std::endl;
     mpf_class *randAGPU = new mpf_class[nr_rows_A * nr_cols_A];
     mpf_class *randAGPUCopy = new mpf_class[nr_rows_A * nr_cols_A];
     std::cout << "Allocated..." << std::endl;
     generateRandomGMPMatrix(randAGPU, nr_rows_A, nr_cols_A);
     int sizeOfArray;
     int maxExp;
     int ownLimbSize = 10;
     int maxAllocSize  = estimateMaxGPUAllocation(0.7, nr_cols_A, nr_cols_B); 
     std::cout << "Maximum size to be allocated is..." << maxAllocSize << std::endl;
     estimateSize(randAGPU, sizeOfArray, maxExp, nr_rows_A, nr_cols_A, ownLimbSize);
     std::cout << "Size of array..." << sizeOfArray << std::endl;
     //std::cout << "Maximum exponent..." << maxExp << std::endl;
     if (sizeOfArray > maxAllocSize) {
       std::cout << "Array too large to be allocated in GPU memory all at once" << std::endl;
     }
     // Allocate memory for temporary matrices
     //double *tmpA = (double *)malloc(nr_rows_A * nr_cols_A * sizeof(double));
     // Allocate memory on the GPU
     double *d_aS; 
     cudaMalloc(&d_aS, sizeOfArray * nr_rows_A * nr_cols_A * sizeof(double));
     //generateLongMatrixFromGMPMatrix_GPU(randAGPU, d_aS, tmpA, sizeOfArray,
     //					 nr_rows_A, nr_cols_A, maxExp,
     //ownLimbSize);
     // Get result back in order to check it 
     //double *h_aS = (double *)malloc(sizeOfArray * nr_rows_A * nr_cols_A * sizeof(double));;
     //cudaMemcpy(h_aS, d_aS, sizeOfArray * nr_rows_A * nr_cols_A * sizeof(double),cudaMemcpyDeviceToHost);
     //std::cout << randA[0].get_str(exp, 2) << std::endl;
     //std::cout << sizeOfArray << std::endl;
     //std::cout << h_aS[10] << std::endl;
     //for(int i = 0; i < sizeOfArray; i++) {
     //  std::bitset<64> trElement(h_aS[i * nr_cols_A * nr_rows_A]);
     //  std::cout << trElement << std::endl;
     //}
     //longToGMPMatrixDouble(randAGPUCopy, h_aS, sizeOfArray, nr_rows_A, nr_cols_A, ownLimbSize, maxExp); 
     //printGMPMatrixDiff(randAGPUCopy, randAGPU, nr_rows_A, nr_cols_A);
     std::cout << "*****************************" << std::endl;
     std::cout << "*** TESTING GMP MATRIX MULTIPLICATION ON GPU ***" << std::endl;
     cublasHandle_t handle;
     cublasCreate(&handle);
     
     // Allocate the rest of the needed arrays to the GPU memory
     double *d_bS; 
     cudaMalloc(&d_bS, sizeOfArray * nr_rows_B * nr_cols_B * sizeof(double));
     double *d_prodRes; 
     cudaMalloc(&d_prodRes, nr_rows_C * nr_cols_C * sizeof(double));
     long  long *d_res; 
     cudaMalloc(&d_res, nr_rows_C * nr_cols_C * sizeof(long long));

     print_memory();
     // Allocate the memory for temporary arrays
     double * tmpTransferLongToGMP = (double *)malloc(max(nr_rows_A * nr_cols_A, nr_rows_B * nr_cols_B) * sizeof(double));
     int * sizeMatrix = (int *)malloc(max(nr_rows_A * nr_cols_A, nr_rows_B * nr_cols_B) * sizeof(int)); 
     int * realExpMatrix = (int *)malloc(max(nr_rows_A * nr_cols_A, nr_rows_B * nr_cols_B) * sizeof(int));
     int * signMatrix = (int *)malloc(max(nr_rows_A * nr_cols_A, nr_rows_B * nr_cols_B) * sizeof(int));
     long long *tmp = (long long *)malloc(nr_rows_C * nr_cols_C * sizeof(long long));
     gettimeofday(&t1, NULL);
     matrixMultiplicationBasecase_cuBlas(handle, randCGPU, randA, randB,
					 d_aS, d_bS, tmpTransferLongToGMP, tmp, sizeMatrix, realExpMatrix, signMatrix,
					 d_prodRes, d_res,  nr_rows_A, nr_cols_A, nr_cols_B);
     gettimeofday(&t2, NULL);
     double etGPU = (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
     //printf("Transfer CPU to GPU memory = %fms\n", etGPU);

   //printGMPMatrix(randCGPU, nr_rows_C, nr_cols_C);
     //printGMPMatrixDiff(randCGPU, randC, nr_rows_C, nr_cols_C);
     printGMPMatrixDiff(randCGPU, randCNaive, 10, 10);//nr_rows_C, nr_cols_C);
     std::cout << std::endl;
     printf("GPU optimized GMP = %fms\n", etGPU);
     printf("CPU naive GMP = %fms\n", etCPU);
     print_memory();
     std::cout << "*****************************" << std::endl;
     // Allocate 3 arrays on CPU    

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
     
     
     std::cout << "Multiplying matrix with transp..." << std::endl;
     // Starting the timer                                                                                           
     gettimeofday(&t1, NULL);

     gpu_blas_mullWithTransp(handle, &d_A[5], d_C, nr_rows_A - 1, nr_cols_A);
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


int main(int argc, char *argv[]) {
  lucaGiantTest(argc, argv);
  //testAddToMpf();
}
