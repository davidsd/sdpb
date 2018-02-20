//
// Created by Petr Kravchuk on 8/15/17.
//

#include <string>
#include <iostream>
#include "../mpmat.h"
#include "mpmat_tests.h"
#include <gmpxx.h>
#include <bitset>
#include <gmp.h>
#include <time.h>

using namespace std;

bool test_mpmatMultiplyGMPBaseCase() {

    cout << "mpmatMultiplyGMPBaseCase(): ";

    const int precs [] = {100, 200, 500, 1000, 2000};
    const int numprecs = 5;
    const int tests = 100;

    gmp_randclass rr(gmp_randinit_default);
    rr.seed(time(NULL));

    for (int i = 0; i < numprecs; i++) {

        int prec = precs[i];

        for(int j=0; j<tests; j++) {
            mpf_class a = 10 * rr.get_f(prec) - 5;
            mpf_class b = 10 * rr.get_f(prec) - 5;
            mpf_class c = 10 * rr.get_f(prec + 256) - 5;

            mpmatMultiplyGMPBaseCase(c, a, b);

            if (!compare_mpf_bits(c, a * b)) { return false; }
        }
    }

    cout << "PASS" << endl;
    return true;
}

