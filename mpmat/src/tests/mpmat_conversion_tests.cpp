//
// Created by Petr Kravchuk on 8/13/17.
//

#include <string>
#include <iostream>
#include "../mpmat.h"
#include <gmpxx.h>
#include <bitset>
#include <gmp.h>
#include <time.h>

using namespace std;

int mpf_test_data_len = 18;
mp_limb_t mpf_test_data [] = {0x3dce6d173c92c9b2, 0x705cc1f9ee6ada0e, 0xbf2e633afc521071, 0x4e537e5bfef78623, 0x1e511ba179a3b41f, 0x6d22a7088f692926, 0x78cadd5b07a608f9, 0xce6d4e0b7064c6ab, 0x8eae7fbce397242c, 0x7c2abfe0c1af9698, 0x5ace25e857cdc622, 0xb558356ae052abff, 0x98d0d40317232ea8, 0xf0b8cdd13d41bb16, 0x4451445d729c05e7, 0x103521440f392c8a, 0x318e60fdcebc4c9b, 0xad23654a1c68};

int mpmat_test_data_26_len = 45;
mp_limb_t mpmat_test_data_26 [] = {0x2b4, 0x23654a1, 0x31a0c63, 0x260fdce, 0x2f1326c, 0x1035214, 0x103ce4b, 0x8a4451, 0x11175ca, 0x1c05e7f, 0x2e3374, 0x13d41bb, 0x5a6343, 0x1403172, 0xcbaa2d, 0x158356a, 0x3814aaf, 0x3f5ace2, 0x17a15f3, 0x1c6227c, 0xaaff83, 0x1af969, 0x223ab9f, 0x3bce397, 0x90b339, 0x2d4e0b7, 0x1931aa, 0x378cadd, 0x16c1e98, 0x8f96d2, 0xa9c223, 0x3692926, 0x79446e, 0x2179a3b, 0x107d394, 0x37e5bfe, 0x3de188e, 0x3f2e633, 0x2bf1484, 0x71705c, 0x307e7b9, 0x2ada0e3, 0x3739b45, 0x33c92c9, 0x2c80000};

int mpmat_test_data_10_len = 116;
mp_limb_t mpmat_test_data_10 [] = {0x0, 0xa, 0x348, 0x365, 0x128, 0x1c6, 0x20c, 0x18e, 0x183, 0x3dc, 0x3af, 0x4c, 0x26c, 0x103, 0x148, 0x144, 0x3c, 0x392, 0x322, 0x244, 0x145, 0x45, 0x35c, 0x29c, 0x17, 0x27f, 0x2e, 0xcd, 0x344, 0x3d4, 0x6e, 0x316, 0x263, 0x10d, 0x100, 0x317, 0x8c, 0x2ea, 0x22d, 0x158, 0xd5, 0x2ae, 0x14, 0x2ab, 0x3fd, 0x1ac, 0x389, 0x1e8, 0x15f, 0xdc, 0x188, 0x27c, 0xaa, 0x3fe, 0x30, 0x1af, 0x25a, 0x188, 0x3ab, 0x27f, 0x2f3, 0x239, 0x1c9, 0x2c, 0x339, 0x2d4, 0x382, 0x370, 0x193, 0x6a, 0x2de, 0xca, 0x375, 0x1b0, 0x1e9, 0x208, 0x3e5, 0x2d2, 0xa9, 0x308, 0x23d, 0x292, 0x249, 0x21e, 0x144, 0x1ba, 0x5e, 0x1a3, 0x2d0, 0x1f4, 0x394, 0x37e, 0x16f, 0x3ef, 0x1e1, 0x223, 0x2fc, 0x2e6, 0xce, 0x2fc, 0x148, 0x107, 0x5c, 0x5c, 0x307, 0x39e, 0x39a, 0x2da, 0x38, 0x3dc, 0x39b, 0x117, 0xf2, 0x12c, 0x26c, 0x200};

int mantissa_size_test_result_sign = -1;
int mantissa_size_test_result_len = 13;
mp_limb_t mantissa_size_test_result [] = {0x8cc7ec153437cc6d, 0x3a629a6cb6e7f83e, 0x636efa377df4cdda, 0x61b244a779bd6ab4, 0xb113fe6d47db2a23, 0x3067471d83e50c1a, 0x7f86084129dd8f08, 0x37c6e14b2443542b, 0xec4b2a257364104a, 0xac0084e589b6fa32, 0xfbe0e21c7036f884, 0xc8e4a385b4b5494, 0x5a88c05676};
int mantissa_size_test_parts_len = 30;
int mantissa_size_test_parts_exp = -52;
mp_limb_signed_t mantissa_size_test_parts [] = {-0x16a2301692b233, 0x3d4bc3d7b9b8b, 0xe8b06630578d, 0x2895874967da2, 0x1eadd0f24861d9, -0x591be5b7e31d5, -0x419e8b36f7f20, 0x156fa3df6d70eb, 0x148baf98e3b781, 0x25f7d5aea503, 0x377dbb5d6afc4, -0xbfe3b79bc385a, -0x1161e7a7c81fbd, -0x433ac869b6478, -0x13aa8f23134e7d, -0x6509efaedc57c, -0x65503e1a531c4, -0x88bcd8c31f48, 0xedc125e013b10, -0x944d31044cc45, -0x10e3996d1ff772, -0x1ee1d476ed7a37, -0x16e7fe1756cb3c, -0x9501943a72c05, -0xd21ddca2cd1a1, -0x5418757544fee, 0x467b5508e6d38, -0x905235a91bb09, -0x13c8b7159b410c, 0x1e80eefbc83393};

bitset<64> print64bits(const void * a){
    return bitset<64>(*static_cast<const long long*>(a));
}

void print_mpf_bits(const mpf_class a) {
    const mp_limb_t * mp_d = a.get_mpf_t() -> _mp_d;
    int size = a.get_mpf_t()-> _mp_size;
    size = size >= 0 ? size : -size;
    cout << "mpf of " << size << " limbs with exp == " << a.get_mpf_t()->_mp_exp << " : " << endl;
    for (int i = size-1; i>=0; i--) {
        cout << "limb " << i << " : " << bitset<64>(mp_d[i]) << endl;
    }
}

void print_mpmat_double_array(const mpmat_double * array, int len) {
    cout << "mpmat_double array of " << len << " limbs" << endl;
    for (int i = 0; i < len; i++) {
        mp_limb_t tmp = static_cast<mp_limb_t>(array[i]);
        cout << "limb " << i << " : " << print64bits(&tmp) << endl;
    }
}

bool test_mpmatConvertGMPToDouble(){

    cout << "mpmatConvertGMPToDouble() : ";

    mpf_t test_number;
    test_number->_mp_d    = mpf_test_data;
    test_number->_mp_prec = mpf_test_data_len - 1;
    test_number->_mp_size = mpf_test_data_len;
    test_number->_mp_exp  = 0;

    mpf_class test_mpf_class (test_number);

    mpmat_double * double_array = new mpmat_double [mpmat_test_data_10_len];

    // Limb of 26 bits

    mpmatConvertGMPToDouble(
            test_mpf_class,
            double_array,
            mpmat_test_data_26_len,
            26,
            0
    );

    bool flag = true;

    for (int i = 0; i < mpmat_test_data_26_len; i++) {
        if (double_array[i] != mpmat_test_data_26[i]) {
            cout << "FAIL during comparison of " << i+1 << "-th limb at mpmat_limb = 26" << endl;

            print_mpf_bits(test_mpf_class);
            print_mpmat_double_array(double_array, mpmat_test_data_26_len);

            flag = false;
            break;
        }
    }

    // Limb of 10 bits

    mpmatConvertGMPToDouble(
            test_mpf_class,
            double_array,
            mpmat_test_data_10_len,
            10,
            0
    );

    for (int i = 0; i < mpmat_test_data_10_len; i++) {
        if (double_array[i] != mpmat_test_data_10[i]) {
            cout << "FAIL during comparison of " << i+1 << "-th limb at mpmat_limb = 10" << endl;

            print_mpf_bits(test_mpf_class);
            print_mpmat_double_array(double_array, mpmat_test_data_10_len);

            flag = false;
            break;
        }
    }

    delete [] double_array;

    if(flag) {
        cout << "PASS" << endl;
    }

    return flag;
}

bool test_mpmatConvertDoubleToGMP(){

    cout << "mpmatConvertDoubleToGMP() : ";

    mpf_t test_number;
    test_number->_mp_d    = mantissa_size_test_result;
    test_number->_mp_prec = mantissa_size_test_result_len - 1;
    test_number->_mp_size = mantissa_size_test_result_sign * mantissa_size_test_result_len;
    test_number->_mp_exp  = 0;

    mpf_class test_result (test_number);
    mpf_class test_parts  ("0", 15*64);

    // Limb of 26 bits

    mpmat_double * double_array = new mpmat_double [mantissa_size_test_parts_len];
    for (int i = 0; i < mantissa_size_test_parts_len; i++) {
        double_array[i] = static_cast<mpmat_double> (mantissa_size_test_parts[i]);
    }

    mpmatConvertDoubleToGMP(
            test_parts,
            double_array,
            mantissa_size_test_parts_len,
            26,
            mantissa_size_test_parts_exp
    );

    bool flag = true;

    for (int i = 0; i < test_result.get_mpf_t()->_mp_size; i++) {
        if (test_result.get_mpf_t()->_mp_d[i] != test_parts.get_mpf_t()->_mp_d[i]) {
            cout << "FAIL during comparison of " << i+1 << "-th limb" << endl;

            print_mpf_bits(test_result);
            print_mpf_bits(test_parts);

            flag = false;
            break;
        }
    }

    delete [] double_array;

    if(flag) {
        cout << "PASS" << endl;
    }

    return flag;
}



bool test_mpmatScalarConversions() {

    cout << "GMP -> Double -> GMP      : ";

    gmp_randclass rr(gmp_randinit_default);
    rr.seed(time(NULL));

    const int num_precs = 5;
    const int precs [] = {100, 200, 500, 1000, 5000};

    const int max_exp  = 128;

    const int mpmat_limbs [] = {16, 26};
    const int num_mpmat_limbs = 2;

    bool flag = true;

    for(int iprec = 0; iprec < num_precs; iprec++) {
        mpf_class r = rr.get_f(precs[iprec])-mpf_class("0.5",precs[iprec]);
        mpf_class t("0",(r.get_mpf_t()->_mp_prec+10)*mp_bits_per_limb);
        t *= 0;

        for(int ilimb = 0; ilimb < num_mpmat_limbs; ilimb++) {
            int size_base = ((r.get_mpf_t()->_mp_prec + 1) * mp_bits_per_limb + max_exp) / mpmat_limbs[ilimb] + 1;
            mpmat_double * double_array = new mpmat_double [size_base];

            for(int exp = 0; exp <= max_exp; exp++) {
                mpmatConvertGMPToDouble(
                        r,
                        double_array,
                        size_base,
                        mpmat_limbs[ilimb],
                        exp
                );

                mpmatConvertDoubleToGMP(
                        t,
                        double_array,
                        size_base,
                        mpmat_limbs[ilimb],
                        exp
                );

                if (r.get_mpf_t()->_mp_size * r.get_mpf_t()->_mp_size < 0) {
                    cout << "FAIL sign mismatch" << endl;
                    cout << "prec: " << precs[iprec] << " limb: " << mpmat_limbs[ilimb] << " exp: " << exp << " diff: " << r-t << endl;
                    print_mpf_bits(r);
                    print_mpf_bits(t);
                    print_mpmat_double_array(double_array, size_base);

                    delete [] double_array;

                    return false;
                }

                for (int i = 0; i < r.get_mpf_t()->_mp_size; i++) {
                    if (r.get_mpf_t()->_mp_d[r.get_mpf_t()->_mp_size-1-i] != t.get_mpf_t()->_mp_d[t.get_mpf_t()->_mp_size-1-i]) {
                        cout << "FAIL during comparison of " << i+1 << "-th most significant limb" << endl;
                        cout << "prec: " << precs[iprec] << " limb: " << mpmat_limbs[ilimb] << " exp: " << exp << " diff: " << r-t << endl;
                        print_mpf_bits(r);
                        print_mpf_bits(t);
                        print_mpmat_double_array(double_array, size_base);

                        delete [] double_array;

                        return false;
                    }
                }
                //print_mpf_bits(r);
                //print_mpf_bits(t);
            }

            delete [] double_array;
        }


    }

    if(flag) {
        cout << "PASS" << endl;
    }

    return flag;

}

void test_run() {

    test_mpmatConvertGMPToDouble();
    test_mpmatConvertDoubleToGMP();
    test_mpmatScalarConversions();

}