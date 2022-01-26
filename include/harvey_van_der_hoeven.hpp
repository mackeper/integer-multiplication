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

    /**
     *  Internal Descrete Fourier Transform, solved by dividing and conquer.
     *  arr: valarray of complex coeffients.
     *  omegas: Complex point used for transformation. (Performance reason)
     *
     *  Can probabilty be optimizied for memory
     */
    void DFT(std::valarray<complex_type> &arr, std::valarray<complex_type> &omegas);

    // transform to unit disc form
    std::vector<complex_type> to_unit_disc(const std::vector<poly_type> &u);

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

    /**
     * Given two intergers (as string) multiply them
     * using Schönhage-Strassen
     */
    int HVDH(std::string fname1, std::string fname2, std::string fout, parameters params);

    /**
    *  Multiplies this polynomial with another, returns a new polynomial.
    *
    *  cannot be const as I need to erase coeffs_
    */
    HVDHPolynomial* HVDHPolynomialMultiplication(HVDHPolynomial &a, HVDHPolynomial &b);

} // namespace imnln
