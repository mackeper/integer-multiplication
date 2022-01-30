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

int main (int argc, char **argv) {
    parse_args(argc, argv);

    srand(0);

    auto params = imnln::get_parameters(imnln::INTEGER_FILE_1, imnln::INTEGER_FILE_2);

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
