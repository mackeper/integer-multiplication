#include <iostream>
#include <fstream>
#include <string>
#include <inttypes.h>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <ctype.h>
#include <stdio.h>
#include <unistd.h>
#include <chrono>

// GMP
#include <gmp.h>

// Local
#include "util.hpp"
#include "tests.hpp"
#include "parameters.hpp"
#include "polynomial_multiplication.hpp"
#include "harvey_van_der_hoeven.hpp"

bool no_gmp = false;
int parse_args(int argc, char **argv) {
    char c;
    while ((c = (char)getopt(argc, argv, "ng:")) != -1)
    {
        printf("%c", c);
        switch (c) {
            case 'g':
                {
                // int val = imnln::MAX_INPUT_SIZE/64;
                int val = imnln::MAX_INPUT_SIZE;
                if (isdigit(optarg[0])) {
                    val = atoi(optarg);
                }
                printf("Generate integer files:\n");

                imnln::timer_start();
                imnln::generate_integer(imnln::INTEGER_FILE_1, val);
                imnln::timer_stop("gen1:\t");

                imnln::timer_start();
                imnln::generate_integer(imnln::INTEGER_FILE_2, val);
                imnln::timer_stop("gen2:\t");
                }
                break;
            case 'n':
                no_gmp = true;
                break;
            /*case 'c':
                cvalue = optarg;
                break;
            case '?':
                if (optopt == 'c')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return 1;*/
            default:
                abort ();
        }
    }
    return 1;
}

int run_gmp(std::string is1,
        std::string is2,
        std::string fgmp,
        std::string fssa,
        std::string fhvdh) {
    imnln::printd("run_gmp start");

    // Multiply
    int flag;
    mpz_t gmp_i1, gmp_i2, gmp_res;
    mpz_init(gmp_i1);
    mpz_set_ui(gmp_i1, 0);
    mpz_init(gmp_i2);
    mpz_set_ui(gmp_i2, 0);
    mpz_init(gmp_res);
    
    {
        std::string is1 = imnln::read_integer(imnln::INTEGER_FILE_1);
        flag = mpz_set_str(gmp_i1, is1.c_str(), 10);
    }
    assert (flag == 0); 
    {
        std::string is2 = imnln::read_integer(imnln::INTEGER_FILE_2);
        flag = mpz_set_str(gmp_i2, is2.c_str(), 10);
    }
    assert (flag == 0); 

    mpz_mul(gmp_res, gmp_i1, gmp_i2);

    mpz_clear(gmp_i1);
    mpz_clear(gmp_i2);

    char * tmp = mpz_get_str(NULL, 10, gmp_res);
    imnln::write_integer(fgmp, tmp);
    delete tmp;

    printf("SSA:  %s\n", imnln::read_integer(fssa).c_str());
    printf("GMP:  %s\n", imnln::read_integer(fgmp).c_str());
    printf("HVDH: %s\n", imnln::read_integer(fhvdh).c_str());
    assert(imnln::compare_integer(fgmp, fssa));
    //assert(imnln::compare_integer(fgmp, fhvdh));
    imnln::printd("run_gmp end");
    return 0;
}

// Change to bin search
uint64_t unique_power_of_two(uint64_t lb, uint64_t ub) {
    for (uint64_t i = lb; i < ub; i++) {
        double log2i = std::log2(i);
        if (std::ceil(log2i) == std::floor(log2i))
            return i;
    }

    return 0;
}

bool is_prime(uint64_t p) {
    for (uint64_t i = 2; i < p; i++) {
        if (p % i == 0) {
            return false;
        }
    }
    return true;
}

uint64_t find_prime_under(uint64_t ub, uint64_t offset) {
    for (uint64_t i = ub-1; i > 0; i--) {
        if (is_prime(i)) {
            if (offset == 0) {
                return i;
            } else {
                offset--;
            }
        }
    }
    return 0;
}

imnln::parameters get_parameters(std::string fname1, std::string fname2) {
    // Read integers from file
    std::string is1 = imnln::read_integer(imnln::INTEGER_FILE_1);
    std::string is2 = imnln::read_integer(imnln::INTEGER_FILE_2);

    // Select parameters
    uint64_t d = 4;
    uint64_t n = 4*std::max(is1.size(), is2.size()); // size of integers
    uint64_t b = (uint64_t)std::ceil(std::log2(n)); // chuck size
    uint64_t p = 6*b; // precision
    // uint64_t a = (uint64_t)std::ceil(std::pow(12*d*d*b, 1/4)); // alpha
    uint64_t a = (uint64_t)std::ceil(std::pow(12*d*d*b, (double)1/4)); // alpha
    uint64_t g = 2*d*a*a; // gamma
    uint64_t T = unique_power_of_two(4*n/b, 8*n/b);
    uint64_t r = unique_power_of_two(
            (uint64_t)std::pow(T, (double)1/(double)d),
            2*(uint64_t)std::pow(T, (double)1/(double)d));
    uint64_t dp = (uint64_t)std::log2(std::pow(r, d)/(double)T);
    std::vector<uint64_t> t(d, r/2); // t
    std::for_each(t.begin()+dp, t.end(), [](uint64_t &i) {i <<= 1;});
    std::vector<uint64_t> s(d); // s
    for (uint64_t i = 0; i < d; i++) {
        if (i < dp) {
            s[i] = find_prime_under(t[i], dp-i-1);
        } else {
            s[i] = find_prime_under(t[i], d-dp-(i-dp)-1);
        }
    }
    uint64_t S = 1;
    std::for_each(s.begin(), s.end(), [&S](uint64_t &i) {S *= i;});
    std::vector<double> o(d); // omegas
    for (uint64_t i = 0; i < o.size(); i++) {
        o[i] = (double)t[i] / (double)s[i]- 1;
        assert(o[i] >= (double)p/std::pow((double)a, (double)4)); // p.36 in paper
    }

    // Set params struct
    imnln::parameters params;
    params.d = d;
    params.n = n;
    params.p = p;
    params.alpha = a;
    params.gamma = g;
    params.T = T;
    params.S = S;
    params.s = s;
    params.t = t;

    printf("Parameters:\n"
            "\td = %" PRIu64 "\n"
            "\tn = %" PRIu64 "\n"
            "\tb = %" PRIu64 "\n"
            "\tp = %" PRIu64 "\n"
            "\ta = %" PRIu64 "\n"
            "\tg = %" PRIu64 "\n"
            "\tT = %" PRIu64 "\n"
            "\tr = %" PRIu64 "\n"
            "\tdp = %" PRIu64 "\n"
            "\tS = %" PRIu64 "\n"
            , d, n, b, p, a, g, T, r, dp, S);
    printf("\n");
    std::for_each(t.begin(), t.end(), [](uint64_t &i) {printf("\tt = %" PRIu64"\n", i);});
    printf("\n");
    std::for_each(s.begin(), s.end(), [](uint64_t &i) {printf("\ts = %" PRIu64"\n", i);});
    printf("\n");
    std::for_each(o.begin(), o.end(), [](double &i) {printf("\to = %lf\n", i);});

    is1.erase();
    is2.erase();

    return params;
}

int main (int argc, char **argv) {
    // imnln::tests();
    parse_args(argc, argv);

    srand(0);

    auto params = get_parameters(imnln::INTEGER_FILE_1, imnln::INTEGER_FILE_2);

    printf("Run multiplications:\n");

    // Run SSA-ish
     imnln::timer_start();
     imnln::SSA(imnln::INTEGER_FILE_1, imnln::INTEGER_FILE_2, imnln::SSA_OUT_FILE);
     imnln::timer_stop("ssa:\t");
    
    imnln::print_ram_info();

    // Run harvey_van_der_hoeven
    imnln::timer_start();
    imnln::HVDH(imnln::INTEGER_FILE_1, imnln::INTEGER_FILE_2, imnln::HVDH_OUT_FILE, params);
    imnln::timer_stop("hvdh:\t");

    imnln::print_ram_info();

    // Compare to gmp (this need lots of ram!)
     if (!no_gmp) {
         imnln::timer_start();
         run_gmp(imnln::INTEGER_FILE_1,
                 imnln::INTEGER_FILE_2,
                 imnln::GMP_OUT_FILE,
                 imnln::SSA_OUT_FILE,
                 imnln::HVDH_OUT_FILE);
         imnln::timer_stop("gmp:\t");
     }

    return 0;
}
