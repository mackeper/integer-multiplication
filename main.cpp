// STD
#include <iostream>
#include <string>
#include <inttypes.h>
#include <cassert>
#include <climits>
#include <cstdlib>

// GMP
#include <gmp.h>

// Local
#include "util.hpp"
#include "tests.hpp"
#include "parameters.hpp"
#include "polynomial_multiplication.hpp"


int main (void) {
    // imnln::tests();

    // std::string is1 = "12032043534543342";
    // std::string is2 = "12032034534534342";
    // std::string is1 = "65536";
    // std::string is2 = "65536";
    // std::string is1 = "4294967296";
    // std::string is2 = "4294967296";
    // std::string is1 = "18446744073709551616";
    // std::string is2 = "18446744073709551616";
    // std::string is1 = "340282366920938463463374607431768211456";
    // std::string is2 = "340282366920938463463374607431768211456";
    // std::string is1 = "115792089237316195423570985008687907853269984665640564039457584007913129639936";
    std::string is2 = "115792089237316195423570985008687907853269984665640564039457584007913129639936";

    std::string is1;
    for (int i = 0; i < 100; i++) {
        is1 += std::to_string(std::rand()%10);       
    }

    std::string ssa_res = imnln::SSA(is1, is2);
    std::cout << ssa_res << std::endl;

    int flag;
    mpz_t gmp_i1, gmp_i2, gmp_res;
    mpz_init(gmp_i1);
    mpz_set_ui(gmp_i1, 0);
    mpz_init(gmp_i2);
    mpz_set_ui(gmp_i2, 0);
    mpz_init(gmp_res);

    flag = mpz_set_str(gmp_i1, is1.c_str(), 10);
    assert (flag == 0); 
    flag = mpz_set_str(gmp_i2, is2.c_str(), 10);
    assert (flag == 0); 

    mpz_mul(gmp_res, gmp_i1, gmp_i2);
    mpz_clear(gmp_i1);
    mpz_clear(gmp_i2);

    char * tmp = mpz_get_str(NULL, 10, gmp_res);
    std::string gmp_res_str = tmp;
    delete tmp;
    mpz_clear(gmp_res);
    std::cout << gmp_res_str << std::endl;

    assert(gmp_res_str.compare(ssa_res) == 0);

    return 0;
}
