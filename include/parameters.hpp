#ifndef INTEGER_MULTIPLICATION_CONSTANTS
#define INTEGER_MULTIPLICATION_CONSTANTS

#include <cstdint>
#include <vector>
#include <string>

namespace imnln {
    const uint64_t DIM = 4;
    const uint64_t MAX_INPUT_SIZE = 1073741824/32;
    const uint64_t CHUCK_SIZE = 4;
    const double EPSILON = CHUCK_SIZE;
    const std::string INTEGER_FILE_1 = "int1.txt";
    const std::string INTEGER_FILE_2 = "int2.txt";
    const std::string SSA_OUT_FILE = "ssa_result.txt";
    const std::string HVDH_OUT_FILE = "hvdh_result.txt";
    const std::string GMP_OUT_FILE = "gmp_result.txt";

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
    parameters get_parameters(const std::string fname1, const std::string fname2);
    uint64_t unique_power_of_two(const uint64_t lb, const uint64_t ub);
    bool is_prime(const uint64_t p);
    uint64_t find_prime_under(const uint64_t ub, uint64_t offset);
}

#endif // INTEGER_MULTIPLICATION_CONSTANTS