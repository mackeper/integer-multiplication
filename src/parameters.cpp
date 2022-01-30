#include <inttypes.h>
#include <algorithm>
#include <cassert>
#include "util.hpp"
#include "parameters.hpp"

// Change to bin search
uint64_t imnln::unique_power_of_two(const uint64_t lb, const uint64_t ub) {
    for (uint64_t i = lb; i < ub; i++) {
        double log2i = std::log2(i);
        if (std::ceil(log2i) == std::floor(log2i))
            return i;
    }
    return 0;
}

bool imnln::is_prime(const uint64_t p) {
    if (p < 2) return false;
    for (uint64_t i = 2; i < p; i++) {
        if (p % i == 0) {
            return false;
        }
    }
    return true;
}

uint64_t imnln::find_prime_under(const uint64_t ub, const uint64_t offset) {
    auto mutable_offset = offset;
    for (uint64_t i = ub-1; i > 0; i--) {
        if (is_prime(i)) {
            if (mutable_offset == 0) {
                return i;
            } else {
                mutable_offset--;
            }
        }
    }
    return 0;
}

imnln::parameters imnln::get_parameters(const std::string fname1, const std::string fname2) {
    // Read integers from file
    //std::string is1 = imnln::read_integer(fname1);
    //std::string is2 = imnln::read_integer(fname2);
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