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
void generateLongVectorFromGMPVector_GPU(const mpf_class *a, long *d_aS, long *tmpA,
                                         int *sizeVector, int *realExpVector, int *signVector, int &sizeOfArray,
                                         const int nr_el, int &maxExp, const int ownLimbSize) {
  assert(ownLimbSize <= INT64L);
  #pragma omp parallel for schedule(dynamic)
  for(int i = 0; i < nr_el; ++i) {
    sizeVector[i] = abs(a[i].get_mpf_t()->_mp_size);
    realExpMatrix[i] = a[i].get_mpf_t()->_mp_exp * INT64L;
    signMatrix[i] =  getSign(a[i]);
  }


  for (int k = 0; k < sizeOfArray - 1; k++) {
    gettimeofday(&t1, NULL);
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < nr_el; ++i) {
      int size = sizeVector[i];
      int realExp = realExpVector[i];
      int sign  = signVector[i];
      //gettimeofday(&tB, NULL);                                                                                 
      //etAccs += (((tB.tv_sec*uS_PER_SEC)+tB.tv_usec) - ((tA.tv_sec*uS_PER_SEC)+tA.tv_usec))/(float)uS_PER_mS; 
      //gettimeofday(&tA, NULL);                                                                                  
 
      int padding = (maxExp - realExp) / ownLimbSize;
      int padBitOffset = (maxExp - realExp) % ownLimbSize;
      //gettimeofday(&tB, NULL);                                                                                  
      //etComputingSizes += (((tB.tv_sec*uS_PER_SEC)+tB.tv_usec) - ((tA.tv_sec*uS_PER_SEC)+tA.tv_usec))/(float)uSER_mSa;                                                                                                               

      if(k > padding) {
	int leftGmpLimb  = size - 1 - (((k - padding) * ownLimbSize - padBitOffset) / INT64L);
	int leftBit      = ((k - padding) * ownLimbSize - padBitOffset) % INT64L;
	int rightGmpLimb = size - 1 - (((k - padding + 1) * ownLimbSize - padBitOffset) / INT64L);
	int rightBit     = ((k - padding + 1) * ownLimbSize - padBitOffset) % INT64L;
	
	if (leftGmpLimb == rightGmpLimb) {
	  long long tmp  = a[i].get_mpf_t()->_mp_d[leftGmpLimb];
	  tmpA[k * nr_el + i] = sign * getBitsFromOneLong(tmp, leftBit, rightBit);
	} else {
	  long long tmp1 = a[i].get_mpf_t()->_mp_d[leftGmpLimb];
	  long long tmp2 = a[i].get_mpf_t()->_mp_d[rightGmpLimb];
	  tmpA[k * nr_el + i] = sign * getBitsFromTwoLong(tmp1, tmp2, leftBit, rightBit);
	}
      } else if (k < padding){
	tmpA[k * nr_el + i] = 0;
      } else {
	long long tmp  = a[i].get_mpf_t()->_mp_d[size - 1];
	tmpA[k * nr_el + i] = sign * getBitsFromOneLong(tmp, 0, ownLimbSize - padBitOffset);
      }
    }
    gettimeofday(&t2, NULL);
    etAlgebra +=  (((t2.tv_sec*uS_PER_SEC)+t2.tv_usec) - ((t1.tv_sec*uS_PER_SEC)+t1.tv_usec))/(float)uS_PER_mS;
    // Transfer matrix to memeory                                                                                   
    //gettimeofday(&t1, NULL);                                                                                      
  }

   #pragma omp parallel for schedule(dynamic)
   for(int i = 0; i < nr_el; ++i) {
     int realExp = realExpVector[i];
     int sign  = signVector[i];

     int padding = (maxExp - realExp) / ownLimbSize;
     int padBitOffset = (maxExp - realExp) % ownLimbSize;
     int leftBit = ((sizeOfArray - padding - 1) * ownLimbSize + padBitOffset) % INT64L;
     int tmp = a[i].get_mpf_t()->_mp_d[0];
     // Make sure to shift this in the appropriate way
     // This is not correct right now but it has a very little effect on the result (only affects the last digits)
     tmpA[(sizeOfArray - 1) * nr_el + i] = sign * getBitsFromOneLong(tmp, leftBit, INT64L);
   }
  cudaMemcpy(d_aS, tmpA, nr_el * sizeOfArray * sizeof(long), cudaMemcpyHostToDevice);
  cudaThreadSynchronize();
  
  //printf("Algebra time in computing hands = %fms\n", etAlgebra);
  //printf("Calculating sizes = %fms\n", etComputingSizes);                                                         
  //printf("My function for obtaining the hands = %fms\n", etAddingLongs);                                          
  //printf("Trivial stuff = %fms\n", etAddTrivial);                                                                 
  //printf("Access time = %fms\n", etAccs);                                                                         
}


// CUDA kernel. Each thread takes care of one element of c                                     
__global__ void vecMult(long *a, long long *b, 
			int nr_el, int sizeOfArray, int nr_rows, long long size)
{
  // Get our global thread ID                                                                                        
  long long id = blockIdx.x*blockDim.x+threadIdx.x;
    
  // Make sure we do not go out of bounds                                                                           
  if (id < size) {
    long id1 = id / (sizeOfArray * nr_el); 
    long id2 = id % (sizeOfArray * nr_el);
    int i = (id1 % sizeOfArray);
    int j = (id2 % sizeOfArray);
    int prec1 = id1/sizeOfArray;
    int prec2 = id2/sizeOfArray;
    if (prec1 + prec2 < sizeOfArray)
      // Replace this with simpler form
      b[sizeOfArray * ((nr_el * (nr_el - 1) - (nr_el - i)  * (nr_el - i -1)) / 2 + j) + prec1 + prec2] += 
	((long long) a[id1]) * a[id2];
  }
}

void subtractVector(long long *b, mpf_class *a, const int nr_el, const int sizeOfArray, const int maxExp, const int ownLimbSize) {
  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < nr_el; ++i) {
    for (int j = 0; j < sizeOfArray; ++j) {
      addToMpf(a[i].get_mpf_t(), b[i * sizeOfArray + j], 2 * maxExp - (j + 2) * ownLimbSize);
    }
  } 
}



void cholesky_GPU(const cublasHandle_t handle, mpf_class *a, 
		  long long *d_aS, long *tmpTransferLongToGMP,
		  long long *tmp, int *sizeMatrix, int *realExpMatrix, int *signMatrix,
		  long long *d_res, long long *d_res,  
		  const int nr_rowsA, const int nr_colsA, const int nr_colsB) {
  // Compute first column of L 
  int sizeOfArray = 0;
  int maxExp = 0;
  int ownLimbSize = INT64L/2 - ceil(log2((double) nr_rows) / 2);
  int tracker = 0;

  estimateSize(a, sizeOfArray, maxExp, ownLimbSize);
  
  mpf_sqrt(a[0], a[0]); 
  
  #pragma omp parallel for schedule(dynamic)
  for (int i = 1; i < nr_rows; ++i) {
    a[i] = a[i]/a[0];
  }
  
  generateLongVectorFromGMPVector_GPU(a, d_aS, tmpTransferLongToGMP,
				      sizeVector, realExpVector, signVector, sizeOfArray,
				      nr_rows  - 1, maxExp, ownLimbSize); 

  // Number of threads in each thread block                                                                                                                 
  int blockSize = 256; // Make sure this number is good? I don't know how                                                                                    

  for(int j = 1; j < nr_cols - GPU_THRESHOLD; ++j) {
    // Number of thread blocks in grid                                            
    int size = sizeOfArray * nr_el;
    size *= size;
    long long gridSize = (long long)ceil((float)(size)/blockSize);
    // Compute results of GPU
    vecMult<<<gridSize, blockSize>>>(d_aS, &d_res[tracker * sizeOfArray],  
				     nr_rows - j, sizeOfArray, nr_rows, size - 1);
    // Copy results from GPU to host
    tracker += nr_rows - j;
    cudaMemcpy(tmp, &d_res[tracker * sizeOfArray], 
	       (nr_rows - j) * sizeOfArray * sizeof(long long), cudaMemcpyDeviceToHost);
    subtractVector(tmp, &a[nr_rows * j + j], nr_rows - j, sizeOfArray);
    mpf_sqrt(a[nr_rows * j + j], a[nr_rows * j + j]);
    
    #pragma omp parallel for schedule(dynamic)
    for(int i = j + 1; i < nr_rows; ++i) {
      a[j * nr_rows + i] = a[j * nr_rows + i]/a[nr_rows * j + j];
    }
    if (j != nr_cols - GPU_THRESHOLD - 1)
      generateLongVectorFromGMPVector_GPU(&a[j * nr_rows + j + 1], d_aS, tmpTransferLongToGMP,
					  sizeVector, realExpVector, signVector, sizeOfArray,
					  nr_rows  - i - 1, maxExp, ownLimbSize);
  }

  // Copy the rest of the results
  tracker += GPU_THRESHOLD;
  cudaMemcpy(tmp, &d_res[tracker * sizeOfArray], (nr_rows * (nr_rows - 1) - tracker) * sizeOfArray * sizeof(long long), cudaMemcpyDeviceToHost);
  for (int j = nr_cols - GPU_THRESHOLD; j < nr_cols; ++j) {
    subtractVector(tmp, &a[j * nr_rows + j], nr_rows - j, sizeOfArray);
    mpf_sqrt(a[j * nr_rows + j], a[j * nr_rows + j]);
    long long diag = a[j * nr_rows + j];
    #pragma omp parallel for schedule(dynamic)
    for (int i = j + 1; i < nr_rows; i++) {
      a[j * nr_rows + i] = a[j * nr_rows + i]/diag; 
      for (int k = j + 1; k <= i; k++) {
	a[k * nr_rows + i] -= a[j * nr_rows + k] * a[j * nr_rows + i];
      }
    }
  }
}
