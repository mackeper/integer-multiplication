#ifndef INTEGER_MULTIPLICATION_UTIL
#define INTEGER_MULTIPLICATION_UTIL
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include "types.hpp"

namespace imnln {
    void timer_start();
    void timer_stop(const std::string text);
    void printd(const std::string s);
    void generate_integer(const std::string fname, const size_t length, const int seed = 0);

    std::string read_integer(const std::string fname);
    bool compare_integer(const std::string fname1, std::string fname2);

    int write_integer(const std::string fname, const poly_type* value, const size_t len, const size_t chuck_size);
    int write_integer(const std::string fname, const char * value);
    poly_type str_to_num(const std::string &istr);
    std::vector<poly_type> split(const std::string &istr, const size_t chuck_size);
    std::string concat(const std::vector<poly_type> &v);

    // source: https://gist.github.com/thirdwing/da4621eb163a886a03c5
    void print_ram_info();

    // Non-square matrix transpose of matrix of size r x c and base address A 
    // https://www.geeksforgeeks.org/inplace-m-x-n-size-matrix-transpose/
    void transpose(std::vector<imnln::poly_type> &v, const int r, const int c);

    // https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/
    // Inverse of a modulo m
    uint64_t mod_inverse(const uint64_t a, const uint64_t m);
}

#endif

