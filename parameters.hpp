#ifndef INTEGER_MULTIPLICATION_CONSTANTS
#define INTEGER_MULTIPLICATION_CONSTANTS

#include <cstdint>
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
}

#endif // INTEGER_MULTIPLICATION_CONSTANTS

