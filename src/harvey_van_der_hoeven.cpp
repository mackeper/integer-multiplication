#include <cmath>
#include <cassert>
#include <valarray>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <fstream>
#include <tuple>
#include <set>
#include <complex>
#include <iomanip>
#include <iostream>
#include "constants.hpp"
#include "types.hpp"
#include "parameters.hpp"
#include "util.hpp"
#include "agarwal_cooley.hpp"
#include "HVDHPolynomial.hpp"
#include "gaussian_resampling.hpp"
#include "harvey_van_der_hoeven.hpp"

/**
 *  Internal Descrete Fourier Transform, solved by dividing and conquer.
 *  arr: valarray of complex coeffients.
 *  omegas: Complex point used for transformation. (Performance reason)
 *
 *  Can probabilty be optimizied for memory
 */
void imnln::DFT(std::valarray<imnln::complex_type> &arr, std::valarray<imnln::complex_type> &omegas) {
    if(arr.size() <= 1) return;

    std::valarray<imnln::complex_type> even = arr[std::slice(0, arr.size()/2, 2)];
    std::valarray<imnln::complex_type> odd = arr[std::slice(1, arr.size()/2, 2)];
    std::valarray<imnln::complex_type> new_omegas = omegas[std::slice(0, arr.size()/2, 2)];
    DFT(even, new_omegas);
    DFT(odd, new_omegas);
    for(size_t i = 0; i < arr.size()/2; i++) {
        auto tmp = omegas[i]*odd[i];
        arr[i] = even[i]+tmp;
        arr[i+arr.size()/2] = even[i]-tmp;
    }
}

// transform to unit disc form
std::vector<imnln::complex_type> imnln::to_unit_disc(const std::vector<imnln::poly_type> &u) {
    std::vector<imnln::complex_type> dv(u.size());
    for (size_t i = 0; i < dv.size(); i++) {
        // this is to unit circle... not disc D:
        dv[i] = imnln::complex_type(std::cos(u[i]), std::sin(u[i]));
    }
    return dv;
}

/**
 * Given two intergers (as string) multiply them
 * using Schönhage-Strassen
 */
int imnln::HVDH(std::string fname1, std::string fname2, std::string fout, imnln::parameters params) {

    // Get first polynomial
    HVDHPolynomial p1;
    {
        std::vector<imnln::poly_type> v1;
        {
            std::string istr1;
            istr1 = imnln::read_integer(fname1);
            v1 = imnln::split(istr1, imnln::CHUCK_SIZE);
            printd("Degree: " + std::to_string(v1.size()) + " S: " 
                    + std::to_string(params.S));
            v1.resize(params.S);
            imnln::agarwal_cooley(v1, params);
            istr1.erase(istr1.begin(), istr1.end());
        }
        p1.vec_to_poly(v1);
    }
    gaussian_resampling(p1, params.S, params.T, params.alpha, params.p);

    // Get second polynomial
    HVDHPolynomial p2;
    {
        std::vector<imnln::poly_type> v2;
        {
            std::string istr2;
            istr2 = imnln::read_integer(fname2);
            v2 = imnln::split(istr2, imnln::CHUCK_SIZE);
            istr2.erase(istr2.begin(), istr2.end());
        }
        p2.vec_to_poly(v2);
    }

    imnln::print_ram_info();

    auto pr = HVDHPolynomialMultiplication(p1, p2);

    // trim 0 coeffients at beginning
    size_t pr_len = pr->size()-1;
    while(pr_len > 0 && std::abs(std::round((*pr)[pr_len])) < imnln::EPSILON) {
        pr_len--;
    }
    // perform carrying
    uint64_t mod = std::pow(10, imnln::CHUCK_SIZE);
    for (size_t i = pr_len; i > 0 && i < __LONG_MAX__; i--) {
        uint64_t v = (uint64_t)std::round((*pr)[i]);
        (*pr)[i-1] += (imnln::poly_type)(v / mod); // carry
        (*pr)[i] = (imnln::poly_type)(v % mod);    // rest
    }

    // Build result string
    imnln::write_integer(fout, pr->get_coeffs(), pr_len, imnln::CHUCK_SIZE);

    delete pr;
    return 0;
}

imnln::HVDHPolynomial* imnln::HVDHPolynomialMultiplication(HVDHPolynomial &a, HVDHPolynomial &b) {
    size_t size = (size_t)std::pow(2,
            std::ceil(std::log2(std::max(a.size(), b.size()))) + 1);

    const std::complex<poly_type> c = 0;

    // Init our complex valarrays
    imnln::printd("Init first complex polynomial");
    std::valarray<std::complex<poly_type>> p1(c, size);
    for(size_t i = 0; i < a.size(); i++)
        p1[i] = std::complex<poly_type>(a[i], 0);
    imnln::print_ram_info();
    a.erase();
    imnln::print_ram_info();

    imnln::printd("Init first complex polynomial");
    std::valarray<std::complex<poly_type>> p2(c, size);
    for(size_t i = 0; i < b.size(); i++)
        p2[i] = std::complex<poly_type>(b[i], 0);
    imnln::print_ram_info();
    b.erase();
    imnln::print_ram_info();

    // Init roots
    imnln::printd("Init roots");
    std::valarray<std::complex<poly_type>> omegas(size);
    for(size_t i = 0; i < size; i++)
        omegas[i] = std::polar((poly_type)1.0, -2.0 * M_PI * ((poly_type)i/(poly_type)size));

    imnln::print_ram_info();

    // Run DFT
    imnln::printd("DFT1");
    imnln::DFT(p1, omegas);
    imnln::printd("DFT2");
    imnln::DFT(p2, omegas);

    imnln::print_ram_info();

    // Dot multiplication
    imnln::printd("point-wise multiplication");
    for(size_t i = 0; i < size; i++)
        p1[i] *= p2[i];

    // Inverse DFT
    imnln::printd("inverse DFT");
    for(size_t i = 0; i < size; i++)
        omegas[i] = std::polar((poly_type)1.0, 2.0 * M_PI * ((poly_type)i/(poly_type)size));
    DFT(p1, omegas);

    imnln::printd("result polynomial");
    auto p = new HVDHPolynomial(size);
    imnln::print_ram_info();
    for(size_t i = 0; i < size; i++) {
        (*p)[i] = p1[i].real()/(poly_type)size;
    }
    return p;
}