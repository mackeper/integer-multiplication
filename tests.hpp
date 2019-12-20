#ifndef INTEGER_MULTIPLICATION_TESTS
#define INTEGER_MULTIPLICATION_TESTS
#include <cstdint>
#include <string>
#include <vector>
#include <cassert>

namespace imnln {
    void tests() {
        // str_to_int
        std::string is0 = "12032043534543342";
        uint64_t i0 = 12032043534543342;
        assert(imnln::str_to_num<uint64_t>(is0) == i0);
        double d0 = 12032043534543342;
        assert(imnln::str_to_num<double>(is0) == d0);

        // split
        std::string is1 = "12032043534543342";
        auto v = imnln::split<uint64_t>(is1);
        assert(is1 == imnln::concat<uint64_t>(v));
        auto vd = imnln::split<double>(is1);
        assert(is1 == imnln::concat<double>(vd));
    }

}

#endif // INTEGER_MULTIPLICATION_TESTS
