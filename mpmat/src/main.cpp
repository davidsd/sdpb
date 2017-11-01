#include <iostream>
#include <bitset>
#include <gmpxx.h>
#include "Timers.h"
#include "mpmat.h"
#include "tests/mpmat_tests.h"
#include "mblas.h"
#include <mkl.h>
#include <omp.h>

Timers timers;

using namespace std;



void print_mpf_info(){
    cout << "sizeof(mp_limb_t)*8 == " << sizeof(mp_limb_t)*8 << endl;
}


template <typename T>
inline T ceil_div(T a, T b){
    return a / b + ( (a % b != 0) & (a > 0) );
}


int main() {

    cout << endl << "============================" << endl;

    print_mpf_info();


    cout << "========== Tests ===========" << endl;

    test_run();

    cout << "============================" << endl;

    mkl_set_num_threads(16);
    omp_set_num_threads(16);

    for (int dim = 2000; dim <= 3000; dim += 1000) {

        cout << " >>>>> Doing dimension " << dim << "<<<<<" << endl;

        int prec = 1000;
        mpf_set_default_prec(prec);

        mpf_class *mat_a = randomGMPVector(dim * dim, prec);
        mpf_class *mat_b = randomGMPVector(dim * dim, prec);
        mpf_class *mat_c = randomGMPVector(dim * dim, prec+256);
        mpf_class *mat_c2 = randomGMPVector(dim * dim, prec);
        mpf_class alpha("1",prec);
        mpf_class beta("0",prec);

        mpmat_gemm_reduced(
                CblasRowMajor, dim, dim, dim,
                mat_a, mat_b,
                mat_c
        );

	mpmat_syrk_reduced(CblasRowMajor, CblasUpper, dim, dim, mat_a, mat_c);

        timers["RgemmParallel"].start();
        RgemmParallel(
                "N", "N", dim, dim, dim,
                alpha,
                mat_a, dim,
                mat_b, dim,
                beta,
                mat_c2,
                dim
        );
        timers["RgemmParallel"].stop();

        cout << timers;

//        for (int i = 0; i < dim*dim; i++) {
//            cout << mat_c[i] - mat_c2[i] << endl;
//        }


//        cout << "Analysis of first element" << endl;
//        cout << mat_c[0]-mat_c2[0] << endl;
//        print_mpf_bits(mat_c[0]);
//        print_mpf_bits(mat_c2[0]);

        delete[] mat_a;
        delete[] mat_b;
        delete[] mat_c;
        delete[] mat_c2;
    }
}
