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
#include "gaussian_resampling.hpp"

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