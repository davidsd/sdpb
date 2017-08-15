#include <iostream>
#include <bitset>
#include <gmpxx.h>
#include "Timers.h"
#include "mpmat.h"
#include "tests/mpmat_tests.h"
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

    mpf_class a("3.1",1024);
    mpf_class b("2.3",1024);
    mpf_class c("0",1024);

    mpmatMultiplyGMPBaseCase(c,a,b);

    cout << "a  : " << a << endl;
    cout << "b  : " << b << endl;
    cout << "a*b: " << c << endl;

    cout << mkl_get_max_threads() << endl;
    cout << omp_get_max_threads() << endl;

    print_mpf_bits(a*b);
    print_mpf_bits(c);
}