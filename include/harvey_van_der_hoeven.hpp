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
