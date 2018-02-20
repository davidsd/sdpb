//
// Created by Petr Kravchuk on 8/13/17.
//

#ifndef MPMAT_MPMAT_TESTS_H
#define MPMAT_MPMAT_TESTS_H


void print_mpf_bits(const mpf_class a);
bool compare_mpf_bits(const mpf_class a, const mpf_class b);

mpf_class * randomGMPVector(int size, int prec);

void print_mpmat_double_array(const mpmat_double * array, int len);

bool test_mpmatConvertGMPToDouble();
bool test_mpmatConvertDoubleToGMP();

// These functionis test invertibility of conversions
bool test_mpmatScalarConversions();
bool test_mpmatVectorConversions();

bool test_mpmatMultiplyGMPBaseCase();

// This function runs all tests
void test_run();

#endif //MPMAT_MPMAT_TESTS_H
