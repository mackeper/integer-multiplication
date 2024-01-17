#ifndef INTEGER_MULTIPLICATION_PARAMETERS
#define INTEGER_MULTIPLICATION_PARAMETERS

#include <cstdint>
#include <vector>
#include <string>

namespace imnln {
    struct parameters {
        uint64_t d; // dimension
        uint64_t n; // size of integers
        uint64_t b; // chuck size
        uint64_t p; // precision
        uint64_t alpha;
        uint64_t gamma;
        uint64_t T;
        uint64_t r;
        uint64_t S;
        std::vector<uint64_t> t;
        std::vector<uint64_t> s;
    };
    parameters get_parameters(const std::string fname1, const std::string fname2, const bool output = false);
    uint64_t unique_power_of_two(const uint64_t lb, const uint64_t ub);
    bool is_prime(const uint64_t p);
    uint64_t find_prime_under(const uint64_t ub, uint64_t offset);
}

#endif