//
// Created by Petr Kravchuk on 8/13/17.
//

#ifndef MPMAT_MPMAT_TESTS_H
#define MPMAT_MPMAT_TESTS_H


void print_mpf_bits(const mpf_class a);

void print_mpmat_double_array(const mpmat_double * array, int len);

bool test_mpmatConvertGMPToDouble();
bool test_mpmatConvertDoubleToGMP();

// This function test invertibility of conversions
bool test_mpmatScalarConversions();

// This function runs all tests
void test_run();

#endif //MPMAT_MPMAT_TESTS_H
