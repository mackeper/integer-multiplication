#ifndef INTEGER_MULTIPLICATION_UTIL
#define INTEGER_MULTIPLICATION_UTIL
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include "types.hpp"

#define DEBUG 1

namespace imnln {
    void timer_start();
    void timer_stop(std::string text);
    void printd(std::string s);
    void generate_integer(std::string fname, size_t length);

    std::string read_integer(std::string fname);
    bool compare_integer(std::string fname1, std::string fname2);

    int write_integer(std::string fname,  poly_type* value, size_t len);
    int write_integer(std::string fname,  char * value);
    poly_type str_to_num(std::string &istr);
    std::vector<poly_type> split(std::string &istr);
    std::string concat(std::vector<poly_type> &v);

    // source: https://gist.github.com/thirdwing/da4621eb163a886a03c5
    void print_ram_info();

    // Non-square matrix transpose of matrix of size r x c and base address A 
    // https://www.geeksforgeeks.org/inplace-m-x-n-size-matrix-transpose/
    template <class T>
    void transpose(std::vector<T> &v, int r, int c);

    // https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/
    // Inverse of a modulo m
    uint64_t modInverse(uint64_t a, uint64_t m);
}

#endif // INTEGER_MULTIPLICATION_UTIL

