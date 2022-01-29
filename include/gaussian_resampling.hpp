#pragma once
#include <cmath>
#include <valarray>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <fstream>
#include <tuple>
#include <set>
#include <complex>
#include <iomanip>
#include "types.hpp"
#include "HVDHPolynomial.hpp"

namespace imnln {
    // S tranformation C^s -> C^t
    std::vector<complex_type> gaussian_resampling_S(const std::vector<complex_type> &u,
            size_t s, size_t t, size_t a, size_t p);

    // Pt transformation, permutation
    std::vector<complex_type> gaussian_resampling_Pt(const std::vector<complex_type> &u,
            size_t s, size_t t);

    // C transformation, "skipping t-s unwanted entries"
    std::vector<complex_type> gaussian_resampling_C(const std::vector<complex_type> &u,
            size_t s, size_t t);

    // I transformation
    std::vector<complex_type> gaussian_resampling_I(const std::vector<complex_type> &u,
            size_t s, size_t t, size_t a, size_t p);

    // D transformation
    std::vector<complex_type> gaussian_resampling_D(const std::vector<complex_type> &u,
            size_t s, size_t t, size_t a);

    // Ps inverse transformation
    std::vector<complex_type> gaussian_resampling_Psinv(const std::vector<complex_type> &u,
            size_t s, size_t t);

    std::vector<complex_type> gaussian_resampling_B(const std::vector<complex_type> &u,
            size_t s, size_t t, size_t a, size_t p);

    int gaussian_resampling(HVDHPolynomial &pv,
            size_t s, size_t t, size_t a, size_t p);
} // namespace imnln
