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
#include "types.hpp"
#include "parameters.hpp"
#include "util.hpp"
#include "agarwal_cooley.hpp"
#include "HVDHPolynomial.hpp"
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

// S tranformation C^s -> C^t
std::vector<imnln::complex_type> imnln::gaussian_resampling_S(const std::vector<imnln::complex_type> &u,
        size_t s, size_t t, size_t a, size_t p) {
    size_t m = (size_t)std::ceil(std::sqrt((imnln::poly_type)p)) * a; 
    m = 10;
    std::vector<imnln::complex_type> tv(t, 0);
    std::vector<std::vector<imnln::poly_type>> mat(t, std::vector<poly_type>(s, 0.0));

    for (size_t k = 0; k < t; k++) {
        // |j - (sk)/t| < m -> (-m + s*k/t, m + s*k/t) 
        int jstart = (int)std::ceil(-(imnln::poly_type)m+((imnln::poly_type)s*(imnln::poly_type)k)/(imnln::poly_type)t);
        int jend   = (int)std::floor((imnln::poly_type)m+((imnln::poly_type)s*(imnln::poly_type)k)/(imnln::poly_type)t);
        for (int j = jstart; j <= jend; j++) {
            imnln::complex_type b = ((imnln::poly_type)1/(imnln::poly_type)2) * ((imnln::poly_type)1/(imnln::poly_type)a);
            //imnln::complex_type b = ((imnln::poly_type)1/(imnln::poly_type)a);
            imnln::complex_type x = -M_PI * std::pow((imnln::poly_type)a, -2);
            x *= std::pow((imnln::poly_type)(j - ((imnln::poly_type)s * (imnln::poly_type)k)/(imnln::poly_type)t), 2);
            x = std::exp(x);

            imnln::complex_type y =  b * x;
            mat[k][((size_t)(j+s*m)) % s] += std::real(y);

            imnln::complex_type z = y * u[((size_t)(j+s*m)) % s];
            tv[k] += z;
        }
    }

    std::cout << std::scientific;
    std::cout << std::setprecision(4);
    for (auto outer : mat) {
        for (auto inner : outer)
            std::cout << inner << " ";
        std::cout << std::endl;
    }
    return tv;
}

// Pt transformation, permutation
std::vector<imnln::complex_type> imnln::gaussian_resampling_Pt(const std::vector<imnln::complex_type> &u,
        size_t s, size_t t) {
    std::vector<imnln::complex_type> ptv(t, 0);
    for (size_t k = 0; k < t; k++) {
        // +(imnln::poly_type)t*(imnln::poly_type)t), because s < t, k < t t*t will make it positive
        // ptv[k] = u[(size_t)((-(imnln::poly_type)s*(imnln::poly_type)k)+(imnln::poly_type)t*(imnln::poly_type)t) % t];
        
        // -sk mod t = (t-s)k mod t?
        ptv[k] = u[(t-s)*k % t];
    }

    return ptv;
}

// C transformation, "skipping t-s unwanted entries"
std::vector<imnln::complex_type> imnln::gaussian_resampling_C(const std::vector<imnln::complex_type> &u,
        size_t s, size_t t) {

    std::vector<imnln::complex_type> csv(s, 0);
    for (size_t l = 0; l < s; l++) {
        csv[l] = u[(size_t)(std::round((imnln::poly_type)t*(imnln::poly_type)l)/(imnln::poly_type)s) % s];
    }
    return csv;
}

// I transformation
std::vector<imnln::complex_type> imnln::gaussian_resampling_I(const std::vector<imnln::complex_type> &u,
        size_t s, size_t t, size_t a, size_t p) {
    std::vector<imnln::complex_type> isv(s, 0);
    std::vector<imnln::complex_type> visv(s, 0);

    //auto ro = [](imnln::poly_type x) -> imnln::poly_type {
    //    return x >= 0 ? std::floor(x) : std::ceil(x);
    //};
    //auto c_ro = [ro](imnln::complex_type x) -> imnln::complex_type {
    //    return imnln::complex_type(ro(std::real(x)), ro(std::imag(x)));
    //};

    // eps transformation
    auto eps = [](const std::vector<imnln::complex_type> &u, size_t s, size_t t, size_t a, size_t p) {
        auto beta = [t, s](size_t l) -> imnln::poly_type {
            return (imnln::poly_type)((imnln::poly_type)t*(imnln::poly_type)l)/(imnln::poly_type)s -
                std::round(((imnln::poly_type)t*(imnln::poly_type)l)/(imnln::poly_type)s);
        };
        
        std::vector<std::vector<poly_type>> mat(s, std::vector<poly_type>(s, 0.0));
        imnln::poly_type m = std::ceil(std::sqrt(p)/((poly_type)(2*a)));
        std::vector<imnln::complex_type> epssv(s, 0);
        for (size_t l = 0; l < s; l++) {
            for (imnln::poly_type h = -(imnln::poly_type)m; h <= (imnln::poly_type)m; h++) {
                if (h == 0) continue; // h not equal 0

                imnln::complex_type x = 1;
                x = -M_PI * std::pow(a, 2);
                x *= (std::pow((imnln::poly_type)t*h/(imnln::poly_type)s + beta(l), 2) -
                        std::pow(beta((size_t)(l+h+s) % s), 2));
                x = std::exp(x);
                mat[l][((size_t)(l+h+s*m)) % s] += std::real(x);

                imnln::complex_type z = x * u[((size_t)(l+h+s*m)) % s];
                epssv[l] += z;
            }
        }
        // std::cout << std::scientific;
        // std::cout << std::setprecision(4);
        // for (auto outer : mat) {
        //     for (auto inner : outer)
        //         std::cout << inner << " ";
        //     std::cout << std::endl;
        // }
        // std::cout << "\n";
        return epssv;
    };

    // v = u/2 page: 30
    for (size_t i = 0; i < s; i++) {
        isv[i] = u[i]/(imnln::poly_type)2;
        visv[i] = u[i]/(imnln::poly_type)2;
        // imnln::complex_type tmp = ((imnln::poly_type)std::pow(2, -10)*c_ro((imnln::poly_type)std::pow(2,9)*u[i]));
        // std::cout << isv[i] << " ro: " << tmp << "\n";
    }

    imnln::poly_type n = std::ceil(p*s/(std::pow(a,2)*(poly_type)(t-s)));
    imnln::poly_type sign = 1;
    for (size_t i = 0; i < n; i++) { // i
        sign *= -1;
        visv = eps(visv, s, t, a, p);
        for (size_t i1 = 0; i1 < s; i1++) { // i
            isv[i1] += sign*visv[i1];
        }
    }

    return isv;
}

// D transformation
std::vector<imnln::complex_type> imnln::gaussian_resampling_D(const std::vector<imnln::complex_type> &u,
        size_t s, size_t t, size_t a) {
    auto beta = [t, s](size_t l) -> imnln::poly_type{
        return (imnln::poly_type)((imnln::poly_type)t*(imnln::poly_type)l)/(imnln::poly_type)s -
            std::round(((imnln::poly_type)t*(imnln::poly_type)l)/(imnln::poly_type)s);
    };

    std::vector<imnln::complex_type> dsv(s, 0);
    for (size_t l = 0; l < s; l++) {
        imnln::complex_type x = 1;
        x = M_PI * std::pow(a, 2) * std::pow(beta(l), 2);
        x = std::exp(x);
        x /= std::pow(2, 2 * std::pow(a, 2) - 2);
        dsv[l] = u[l]*x;
    }

    return dsv;
}

// Ps inverse transformation
std::vector<imnln::complex_type> imnln::gaussian_resampling_Psinv(const std::vector<imnln::complex_type> &u,
        size_t s, size_t t) {
    std::vector<imnln::complex_type> pssv(s, 0);
    std::vector<std::tuple<size_t, imnln::complex_type>> tpssv(s);

    for (size_t i = 0; i < s; i++) {
        tpssv[i] = std::make_pair((t*i) % s, u[i]);
    }

    auto complex_cmp = [](const std::tuple<size_t, imnln::complex_type> &t1,
            const std::tuple<size_t, imnln::complex_type> &t2) {
        return (std::get<0>(t1) == std::get<0>(t2)
            && std::real(std::get<1>(t1)) < std::real(std::get<1>(t2)))
            || std::real(std::get<0>(t1)) < std::real(std::get<0>(t2));
    };

    std::sort(tpssv.begin(), tpssv.end(), complex_cmp);

    for (size_t i = 0; i < s; i++) {
        pssv[i] = std::get<1>(tpssv[i]);
    }

    return pssv;
}

std::vector<imnln::complex_type> imnln::gaussian_resampling_B(const std::vector<imnln::complex_type> &u,
        size_t s, size_t t, size_t a, size_t p) {
    imnln::printd("Start gaussian resampling inv");
    std::vector<imnln::complex_type> ub(u);
    ub = gaussian_resampling_Pt(ub, s, t);
    ub = gaussian_resampling_C(ub, s, t);
    ub = gaussian_resampling_I(ub, s, t, a, p);
    ub = gaussian_resampling_D(ub, s, t, a);
    ub = gaussian_resampling_Psinv(ub, s, t);
    imnln::printd("End gaussian resampling inv");
    return ub;
}

int imnln::gaussian_resampling(HVDHPolynomial &pv,
        size_t s, size_t t, size_t a, size_t p) {
    std::cout << std::setprecision(5);
    std::cout << std::fixed;
    imnln::printd("Start gaussian resampling");

    // imnln::poly_type a = (imnln::poly_type)imnln::params.alpha;
    // imnln::poly_type s = (imnln::poly_type)imnln::params.s.size();
    // imnln::poly_type t = (imnln::poly_type)imnln::params.t.size();
    
    s = 10;
    t = 13;
    a = 2;
    p = 100;
    assert((imnln::poly_type)t/(imnln::poly_type)s - 1.0 > 1.0/((imnln::poly_type)a*(imnln::poly_type) a));

    // std::vector<imnln::poly_type> rsv = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    std::vector<imnln::poly_type> rsv = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    for (auto &i : rsv) {
        std::cout << i << std::endl;
    }
    std::cout << std::endl;
    // std::vector<imnln::complex_type> sv = to_unit_disc(rsv);
    std::vector<imnln::complex_type> sv(rsv.size());
    for (size_t i = 0; i < sv.size(); i++) {
        // sv[i] = 0.01 * rsv[i];
        sv[i] = (imnln::poly_type)std::pow(2, -(imnln::poly_type)p)*rsv[i];
        // sv[i] = rsv[i];
        std::cout << sv[i] << "\n";
    }
    std::cout << std::endl;
    std::vector<imnln::complex_type> tv = gaussian_resampling_S(sv, s, t, a, p);

    for (auto &i : tv) {
        std::cout << i << std::endl;
    }
    std::cout << std::endl;
    imnln::printd("End gaussian resampling");


    // const imnln::complex_type c = 0;
    // std::valarray<imnln::complex_type> p1(c, t);
    // for(size_t i = 0; i < t; i++) {
    //     p1[i] = imnln::complex_type(std::real(tv[i]), 0);
    // }
    // // Init roots
    // std::valarray<imnln::complex_type> omegas(t);
    // for(size_t i = 0; i < t; i++) {
    //     omegas[i] = std::polar((imnln::poly_type)1.0, -2.0 * M_PI * ((imnln::poly_type)i/(imnln::poly_type)t));
    // }
    // DFT(p1, omegas); // Run DFT
    // for(size_t i = 0; i < t; i++) {
    //     omegas[i] = std::polar((imnln::poly_type)1.0, 2.0 * M_PI * ((imnln::poly_type)i/(imnln::poly_type)t));
    // }
    // DFT(p1, omegas); // Inverse DFT
    // std::vector<imnln::complex_type> tvdft(t);
    // for (size_t i = 0; i < t; i++) {
    //     tvdft[i] = p1[i]/(imnln::poly_type)t;
    //     std::cout << tvdft[i] << std::endl;
    // }

    auto pssv = gaussian_resampling_B(tv, s, t, a, p);
    for (size_t i = 0; i < s; i++) {
        std::cout << std::real((pssv[i]*(imnln::complex_type)std::pow(2, 2*a*a)))*2 << "   ";
        std::cout << std::real(sv[i]) << "    ";
        std::cout << std::real((pssv[i]*(imnln::complex_type)std::pow(2, 2*a*a))/sv[i]) << std::endl;
    }
    std::cout << std::endl;
    return 0;
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
            istr1 = imnln::read_integer(imnln::INTEGER_FILE_1);
            v1 = imnln::split(istr1);
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
            istr2 = imnln::read_integer(imnln::INTEGER_FILE_2);
            v2 = imnln::split(istr2);
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
    imnln::write_integer(fout, pr->get_coeffs(), pr_len);

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