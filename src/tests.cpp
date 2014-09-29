#include <vector>
#include <iostream>
#include "types.h"
#include "Vector.h"
#include "Matrix.h"
#include "BlockDiagonalMatrix.h"

using std::vector;
using std::cout;
using std::endl;

void testCholeskyStabilize() {
  cout.precision(20);
  Matrix A(4,4);
  Matrix L(A);
  A.elt(0,0) = 1e20;
  A.elt(1,1) = 2e18;
  A.elt(2,2) = 1e-3;
  A.elt(3,3) = 1e-5;
  A.elt(1,0) = 1e2;
  A.elt(0,1) = 1e2;
  A.elt(1,2) = 2e2;
  A.elt(2,1) = 2e2;
  vector<int> updateIndices;
  Vector updateVector(L.rows);
  Real lambdaGM;
  choleskyDecomposition(A,L);
  cout << "A = " << A << ";\n";
  cout << "L = " << L << ";\n";
  stabilizeCholesky(L, updateVector, updateIndices, lambdaGM);
  cout << "updateIndices = " << updateIndices << ";\n";
  cout << "L = " << L << "\n;";
  cout << "lambdaGM = " << lambdaGM << ";\n";


  vector<Integer> stabilizeIndices;
  vector<Real> stabilizeLambdas;
  choleskyDecompositionStabilized(A, L, stabilizeIndices, stabilizeLambdas);
  cout << "A = " << A << ";\n";
  cout << "L = " << L << ";\n";
  cout << "stabilizeIndices = " << stabilizeIndices << ";\n";
  cout << "stabilizeLambdas = " << stabilizeLambdas << ";\n";
}

void testLinearlyIndependentRowIndices() {
  Matrix A(8,3);
  A.elt(0,0)=0;
  A.elt(0,1)=0;
  A.elt(0,2)=0;
  A.elt(1,0)=1;
  A.elt(1,1)=0;
  A.elt(1,2)=1;
  A.elt(2,0)=0;
  A.elt(2,1)=1;
  A.elt(2,2)=0;
  A.elt(3,0)=0;
  A.elt(3,1)=0;
  A.elt(3,2)=0;
  A.elt(4,0)=0;
  A.elt(4,1)=1;
  A.elt(4,2)=1;
  A.elt(5,0)=1;
  A.elt(5,1)=1;
  A.elt(5,2)=0;
  A.elt(6,0)=1;
  A.elt(6,1)=0;
  A.elt(6,2)=1;
  A.elt(7,0)=1;
  A.elt(7,1)=0;
  A.elt(7,2)=1;

  cout << "A = " << A << ";\n";
  
  vector<int> rows = linearlyIndependentRowIndices(A);
  
  cout << "Aprime = " << A << ";\n";
  cout << "rows = " << rows << ";\n";
}

void testCholeskyUpdate() {
  Matrix A(4,4);
  Matrix B(A);
  Matrix C(A);
  Matrix L(A);
  Matrix LT(L);
  Matrix V(4, 2);
  Matrix VT(V.cols, V.rows);
  V.elt(0,0) =1;
  V.elt(1,0) =2;
  V.elt(2,0) =3;
  V.elt(3,0) =4;
  V.elt(0,1) =5;
  V.elt(1,1) =4;
  V.elt(2,1) =3;
  V.elt(3,1) =2;
  for (int r = 0; r < V.rows; r++)
    for (int c = 0; c < V.cols; c++)
      VT.elt(c, r) = V.elt(r,c);
  Matrix U(V);

  A.addDiagonal(4);
  cout << "A = " << A << endl;
  cout << "V = " << V << endl;
  choleskyDecomposition(A, L);
  transpose(L, LT);

  matrixMultiply(V, VT, B);
  B += A;
  matrixMultiply(L, LT, C);
  C -= B;

  cout << "L L^T - (A + V V^T) = " << C << endl;
}

void testMatrix() {
  Matrix A(3,3);
  A.elt(0,0) = 1;
  A.elt(1,0) = 2;
  A.elt(2,0) = 3;
  A.symmetrize();
  cout << A << endl;
}
