#include "../../src/SDPSolver.h"
#include "../../src/Timers.h"
#include "../../src/parse.h"

// FIXME: This should not be a global
Timers timers;

int main() {

  mp_bitcnt_t precision(2048);
  mpf_set_default_prec(precision);

  Real residue("10e-" + std::to_string(std::floor(precision/(log(10)/log(2)))));
  Real sqrt_residue(sqrt(residue));
  
  const size_t N(1);
  Matrix matrix(N, N), square_gmp(N, N);
  for (size_t i = 0; i < N; ++i)
    for (size_t j = 0; j < N; ++j) {
      Real ri(i + 1);
      matrix.elt(i, j) = ri-sqrt_residue;
    }

  matrixSquareIntoBlock(matrix, square_gmp, 0, 0);

  mpmat workspace(1);
  mpf_set_default_prec(mpf_get_default_prec() + 256);
  Matrix square_mpmat(N, N);
  mpf_set_default_prec(mpf_get_default_prec() - 256);

  workspace.syrk_reduced(CblasRowMajor, CblasTrans, matrix.cols, matrix.rows,
                         matrix.elements.data(), square_mpmat.elements.data());

  // matrixSquareIntoBlockMpmat(workspace, matrix, square_mpmat, 0, 0);

  bool passed (true);
  std::cout.precision(static_cast<int>(precision * 0.31 + 5));
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (abs((square_gmp.elt(i, j) - square_mpmat.elt(i, j)) /
              (square_gmp.elt(i, j) + square_mpmat.elt(i, j))) > residue)
        {
          passed=false;
          std::cout << i << " " << j << "\n\t"
                    << "square_gmp: " << square_gmp.elt(i, j) << "\n\t"
                    << "square_mpmat: " << square_mpmat.elt(i, j) << "\n\t"
                    << "target residue: " << residue << "\n\t"
                    << "measured residue: "
                    << (square_gmp.elt(i, j) - square_mpmat.elt(i, j)) << "\n";
        }
    }
  }
  if(passed)
    { std::cout << "PASS: sqrt(1-epsilon)^2\n"; }
  else
    { std::cout << "FAIL: sqrt(1-epsilon)^2\n"; }
  
  // SDPSolverParameters parameters;
  // parameters.precision = 2048;
  // parameters.resetPrecision();
  // parameters.maxThreads = 1;
  // // parameters.maxThreads=128;

  // omp_set_num_threads(parameters.maxThreads);

  // SDPSolver solver(readBootstrapSDP(std::vector<boost::filesystem::path>(
  //                      // {"../inputs/TTTTtestsmall.xml"})),
  //                      {"test/mpmat/input.xml"})),
  //                  parameters);
  // solver.testMultiplication(1, 1, 5);
}
