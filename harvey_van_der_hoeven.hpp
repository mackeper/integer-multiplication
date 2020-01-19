// Authors: Marcus Östling
#pragma once
#include <cmath>
#include <complex>
#include <valarray>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <fstream>
#include <tuple>
#include <set>
#include <complex>

namespace imnln {

    typedef long double poly_type;
    typedef std::complex<poly_type> complex_type;

    /**
     *  Internal Descrete Fourier Transform, solved by dividing and conquer.
     *  arr: valarray of complex coeffients.
     *  omegas: Complex point used for transformation. (Performance reason)
     *
     *  Can probabilty be optimizied for memory
     */
    void DFT(std::valarray<complex_type> &arr, std::valarray<complex_type> &omegas) {
        if(arr.size() <= 1) return;

        std::valarray<complex_type> even = arr[std::slice(0, arr.size()/2, 2)];
        std::valarray<complex_type> odd = arr[std::slice(1, arr.size()/2, 2)];
        std::valarray<complex_type> new_omegas = omegas[std::slice(0, arr.size()/2, 2)];
        DFT(even, new_omegas);
        DFT(odd, new_omegas);
        for(size_t i = 0; i < arr.size()/2; i++) {
            auto tmp = omegas[i]*odd[i];
            arr[i] = even[i]+tmp;
            arr[i+arr.size()/2] = even[i]-tmp;
        }
    }


    /**
     *  Class that stores the coeffients of a polynomial and using
     *  FFT can be multiplied with another HVDHPolynomial.
     */
    template <class T>
    class HVDHPolynomial {

        T* coeffs_;
        size_t size_;
        bool deleted = true;

        public:
        HVDHPolynomial() {}

        explicit HVDHPolynomial(size_t size) {
            coeffs_ = new T[size];
            size_ = size;
            deleted = false;
        }

        explicit HVDHPolynomial(std::vector<T> &v) {
            coeffs_ = new T[v.size()];
            size_ = v.size();
            for (size_t i = 0; i < v.size(); i++) {
                coeffs_[i] = v[i];
            }
            deleted = false;
        }

        HVDHPolynomial(HVDHPolynomial<T> &&p) {
            this->size_ = p.size_;
            this->coeffs_ = p.coeffs_;
            p.coeffs_ = nullptr;
            deleted = false;
        }

        ~HVDHPolynomial() {
            if (!deleted) {
                delete[] coeffs_;
                size_ = 0;
                deleted = true;
            }
        }

        void erase() {
            if (!deleted) {
                delete[] coeffs_;
                size_ = 0;
                deleted = true;
            }
        }

        void vec_to_poly(std::vector<T> &v) {
            this->erase();
            coeffs_ = new T[v.size()];
            size_ = v.size();
            for (size_t i = 0; i < v.size(); i++) {
                coeffs_[i] = v[i];
            }
            deleted = false;
        }

        size_t size() const {
            return size_;
        }

        T& operator[](size_t i) {
            return coeffs_[i];
        }

        T operator[](size_t i) const {
            return coeffs_[i];
        }

        T* get_coeffs(){
            return coeffs_;
        }

        /**
         *  Multiplies this polynomial with another, returns a new polynomial.
         *
         *  cannot be const as I need to erase coeffs_
         */
        HVDHPolynomial<T> polymul(HVDHPolynomial<T> &other) {
            size_t size = (size_t)std::pow(2,
                    std::ceil(std::log2(std::max(size_, other.size()))) + 1);

            const std::complex<T> c = 0;

            // Init our complex valarrays
            imnln::printd("Init first complex polynomial");
            std::valarray<std::complex<T>> p1(c, size);
            for(size_t i = 0; i < size_; i++)
                p1[i] = std::complex<T>((*this)[i], 0);
            imnln::print_ram_info();
            this->erase();
            imnln::print_ram_info();

            imnln::printd("Init first complex polynomial");
            std::valarray<std::complex<T>> p2(c, size);
            for(size_t i = 0; i < other.size(); i++)
                p2[i] = std::complex<T>(other[i], 0);
            imnln::print_ram_info();
            other.erase();

            imnln::print_ram_info();

            // Init roots
            imnln::printd("Init roots");
            std::valarray<std::complex<T>> omegas(size);
            for(size_t i = 0; i < size; i++)
                omegas[i] = std::polar((T)1.0, -2.0 * M_PI * ((T)i/(T)size));

            imnln::print_ram_info();

            // Run DFT
            imnln::printd("DFT1");
            DFT(p1, omegas);
            imnln::printd("DFT2");
            DFT(p2, omegas);

            imnln::print_ram_info();

            // Dot multiplication
            imnln::printd("point-wise multiplication");
            for(size_t i = 0; i < size; i++)
                p1[i] *= p2[i];

            // Inverse DFT
            imnln::printd("inverse DFT");
            for(size_t i = 0; i < size; i++)
                omegas[i] = std::polar((T)1.0, 2.0 * M_PI * ((T)i/(T)size));
            DFT(p1, omegas);

            imnln::printd("result polynomial");
            HVDHPolynomial p(size);
            imnln::print_ram_info();
            for(size_t i = 0; i < size; i++) {
                p[i] = p1[i].real()/(T)size;
            }
            return p;
        }
    };

    int agarwal_cooley(std::vector<poly_type> &v) {
        std::vector<uint64_t> s = imnln::params.s;
        std::vector<uint64_t> t = imnln::params.t;
        uint64_t S = imnln::params.S;
        uint64_t T = imnln::params.T;

        std::vector<uint64_t> v2 = {0, 1, 2, 3, 4, 4, 6, 7, 8, 9};
        std::vector<uint64_t> invs(s.size());
        for (uint64_t i = 0; i < s.size(); i++) {
            invs[i] = imnln::modInverse(t[i], s[i]);
            std::cout << "invs: " << invs[i] << std::endl;
        }

        std::vector<std::tuple<uint64_t, uint64_t>> v3;
        // size_t N = 10;
        // size_t r1 = 2;
        // size_t r2 = 5;
        // size_t c = 2;


        // uint64_t s1 = modInverse(r2, r1) * r2;
        // uint64_t s2 = modInverse(r1, r2) * r1;
        // std::cout << "s1: " << s1 << " s2: " << s2 << "\n";
        // for (auto i : v2) {
        //     uint64_t i1 = i % r1;
        //     uint64_t i2 = i % r2;
        //     for (uint64_t x = 0; x < s.size(); x++) {
        //
        //     }
        //     uint64_t inv = (i1*s1 + i2*s2) % N;
        //     for (auto i : s) {
        //
        //     }
        //
        //     std::cout << "i: " << i << " i1: " << i1 << " i2: " << i2 <<
        //         " inverse: " << inv << "\n";
        //     v3.push_back(std::make_tuple(i1, i2));
        //     assert(i == inv);
        // }
        //
        // for (uint64_t k2 = 0; k2 < r2; k2++) {
        //     for (uint64_t k1 = 0; k1 < r1; k1++) {
        //
        //     }
        // }

        // return 0;
        //std::vector<int> v2 = {1,2,3,4,5,6};
        // for (size_t i = 0; i < n; i++) {
        //     for (size_t j = 0; j < m; j++) {
        //         std::cout << v2[i*m+j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;
        // imnln::transpose(v2, n, m);
        // for (size_t i = 0; i < m; i++) {
        //     for (size_t j = 0; j < n; j++) {
        //         std::cout << v2[i*n+j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        
        // std::cout << "Y: " << std::endl;
        // for (size_t i = 0; i < n; i++)  {
        //     for (size_t j = 0; j < m; j++)  {
        //         // std::cout << v2[(((i-j*c)%m)*m+j)*c] << std::endl;
        //         std::cout << v2[((i-j*c)%m)*m+j] << std::endl;
        //     }
        // }
        //
        // return 0;
        //
        std::vector<std::tuple<uint64_t, poly_type>> label_and_sort(S);
        std::set<uint64_t> st;
        for (size_t i = 0; i < S; i++) {
            label_and_sort[i] = std::tuple<poly_type, poly_type>((T*i) % S, v[i]);
            st.insert((T*i)%S);
        }

        std::sort(label_and_sort.begin(), label_and_sort.end());

        for (size_t i = 0; i < S; i++) {
            // std::cout << std::get<0>(label_and_sort[i]) << ", "
            //     << std::get<1>(label_and_sort[i]) << std::endl;
            if (invs[0] * i < S) {
                std:: cout << invs[0] * i << std::endl;
                // std::cout << std::get<1>(label_and_sort[invs[0]*i]) << std::endl;
            }
        }
        return 0;
    }

    // transform to unit disc form
    std::vector<complex_type> to_unit_disc(const std::vector<poly_type> &u) {
        std::vector<complex_type> dv(u.size());
        for (size_t i = 0; i < dv.size(); i++) {
            // this is to unit circle... not disc D:
            dv[i] = complex_type(std::cos(u[i]), std::sin(u[i]));
        }
        return dv;
    }

    // S tranformation C^s -> C^t
    std::vector<complex_type> gaussian_resampling_S(const std::vector<complex_type> &u,
            size_t s, size_t t, size_t a, size_t p) {
        size_t m = (size_t)std::ceil(std::sqrt((poly_type)p)) * a; 
        std::vector<complex_type> tv(t, 0);

        for (size_t k = 0; k < t; k++) {
            // |j - (sk)/t| < m -> (-m + s*k/t, m + s*k/t) 
            poly_type jstart = -(poly_type)m+((poly_type)s*(poly_type)k)/(poly_type)t + 1;
            poly_type jend   =  (poly_type)m+((poly_type)s*(poly_type)k)/(poly_type)t;

            for (poly_type j = jstart; j < jend; j++) {
                complex_type b = ((poly_type)1/(poly_type)2) * ((poly_type)1/(poly_type)a);

                complex_type x = -M_PI * std::pow((poly_type)a, -2);
                x *= std::pow((poly_type)j - ((poly_type)s * (poly_type)k)/(poly_type)t, 2);
                x = std::exp(x);

                complex_type y =  b * x;

                complex_type z = y * u[(size_t)(j + s) % s];
                tv[k] += z;
            }
        }
        return tv;
    }

    // Pt transformation, permutation
    std::vector<complex_type> gaussian_resampling_Pt(const std::vector<complex_type> &u,
            size_t s, size_t t) {
        std::vector<complex_type> ptv(t, 0);
        for (size_t k = 0; k < t; k++) {
            // +(poly_type)t*(poly_type)t), because s < t, k < t t*t will make it positive
            // ptv[k] = u[(size_t)((-(poly_type)s*(poly_type)k)+(poly_type)t*(poly_type)t) % t];
            
            // -sk mod t = (t-s)k mod t?
            ptv[k] = u[(t-s)*k % t];
        }

        return ptv;
    }

    // C transformation, "skipping t-s unwanted entries"
    std::vector<complex_type> gaussian_resampling_C(const std::vector<complex_type> &u,
            size_t s, size_t t) {
    
        std::vector<complex_type> csv(s, 0);
        for (size_t l = 0; l < s; l++) {
            csv[l] = u[(size_t)(std::round((poly_type)t*(poly_type)l)/(poly_type)s) % s];
        }
        return csv;
    }

    // I transformation
    std::vector<complex_type> gaussian_resampling_I(const std::vector<complex_type> &u,
            size_t s, size_t t, size_t a, size_t p) {
        std::vector<complex_type> isv(s, 0);
        std::vector<complex_type> visv(s, 0);

        auto ro = [](poly_type x) -> poly_type {
            return x >= 0 ? std::floor(x) : std::ceil(x);
        };
        auto c_ro = [ro](complex_type x) -> complex_type {
            return complex_type(ro(std::real(x)), ro(std::imag(x)));
        };

        // eps transformation
        auto eps = [](const std::vector<complex_type> &u, size_t s, size_t t, size_t a, size_t p) {
            auto beta = [t, s](size_t l) -> poly_type {
                return (poly_type)((poly_type)t*(poly_type)l)/(poly_type)s -
                    std::round(((poly_type)t*(poly_type)l)/(poly_type)s);
            };
            
            poly_type m = std::ceil(std::sqrt(p)/(2*a));
            std::vector<complex_type> epssv(s, 0);
            for (size_t l = 0; l < s; l++) {
                for (poly_type h = -(poly_type)m; h <= (poly_type)m; h++) {
                    if (h == 0) continue; // h not equal 0

                    complex_type x = 1;
                    x = -M_PI * std::pow(a, 2);
                    x *= (std::pow((poly_type)t*h/(poly_type)s + beta(l), 2) -
                            std::pow(beta((size_t)(l+h+s) % s), 2));
                    x = std::exp(x);

                    complex_type z = x * u[(size_t)(l+h+s) % s];
                    epssv[l] += z;
                }
            }
            return epssv;
        };

        // v = u/2 page: 30
        for (size_t i = 0; i < s; i++) {
            isv[i] = u[i]/(poly_type)2;
            visv[i] = u[i]/(poly_type)2;
            // complex_type tmp = ((poly_type)std::pow(2, -10)*c_ro((poly_type)std::pow(2,9)*u[i]));
            // std::cout << isv[i] << " ro: " << tmp << "\n";
        }

        poly_type n = std::ceil(p*s/(std::pow(a,2)*(t-s)));
        poly_type sign = 1;
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
    std::vector<complex_type> gaussian_resampling_D(const std::vector<complex_type> &u,
            size_t s, size_t t, size_t a) {
        auto beta = [t, s](size_t l) -> poly_type{
            return (poly_type)((poly_type)t*(poly_type)l)/(poly_type)s -
                std::round(((poly_type)t*(poly_type)l)/(poly_type)s);
        };

        std::vector<complex_type> dsv(s, 0);
        for (size_t l = 0; l < s; l++) {
            complex_type x = 1;
            x = M_PI * std::pow(a, 2) * std::pow(beta(l), 2);
            x = std::exp(x);
            x /= std::pow(2, 2 * std::pow(a, 2) - 2);
            dsv[l] = u[l]*x;
        }

        return dsv;
    }

    // Ps inverse transformation
    std::vector<complex_type> gaussian_resampling_Psinv(const std::vector<complex_type> &u,
            size_t s, size_t t) {
        std::vector<complex_type> pssv(s, 0);
        std::vector<std::tuple<size_t, complex_type>> tpssv(s);

        for (size_t i = 0; i < s; i++) {
            tpssv[i] = std::make_pair((t*i) % s, u[i]);
        }

        auto complex_cmp = [](const std::tuple<size_t, complex_type> &t1,
                const std::tuple<size_t, complex_type> &t2) {
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

    int gaussian_resampling(HVDHPolynomial<poly_type> &pv) {
        std::cout << std::setprecision(5);
        std::cout << std::fixed;
        imnln::printd("Start gaussian resampling");

        // poly_type a = (poly_type)imnln::params.alpha;
        // poly_type s = (poly_type)imnln::params.s.size();
        // poly_type t = (poly_type)imnln::params.t.size();
        
        size_t s = 3;
        size_t t = 4;
        size_t a = 1;
        size_t p = 1;

        // std::vector<poly_type> rsv = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
        std::vector<poly_type> rsv = {1, 1, 1};

        for (auto &i : rsv) {
            std::cout << i << std::endl;
        }
        std::cout << std::endl;
        // std::vector<complex_type> sv = to_unit_disc(rsv);
        std::vector<complex_type> sv(rsv.size());
        for (size_t i = 0; i < sv.size(); i++) {
            // sv[i] = 0.01 * rsv[i];
            sv[i] = (poly_type)std::pow(2, -(poly_type)p)*rsv[i];
            // sv[i] = rsv[i];
            std::cout << sv[i] << "\n";
        }
        std::cout << std::endl;
        std::vector<complex_type> tv = gaussian_resampling_S(sv, s, t, a, p);

        for (auto &i : tv) {
            std::cout << i << std::endl;
        }
        std::cout << std::endl;
        imnln::printd("End gaussian resampling");


        // const complex_type c = 0;
        // std::valarray<complex_type> p1(c, t);
        // for(size_t i = 0; i < t; i++) {
        //     p1[i] = complex_type(std::real(tv[i]), 0);
        // }
        // // Init roots
        // std::valarray<complex_type> omegas(t);
        // for(size_t i = 0; i < t; i++) {
        //     omegas[i] = std::polar((poly_type)1.0, -2.0 * M_PI * ((poly_type)i/(poly_type)t));
        // }
        // DFT(p1, omegas); // Run DFT
        // for(size_t i = 0; i < t; i++) {
        //     omegas[i] = std::polar((poly_type)1.0, 2.0 * M_PI * ((poly_type)i/(poly_type)t));
        // }
        // DFT(p1, omegas); // Inverse DFT
        // std::vector<complex_type> tvdft(t);
        // for (size_t i = 0; i < t; i++) {
        //     tvdft[i] = p1[i]/(poly_type)t;
        //     std::cout << tvdft[i] << std::endl;
        // }

        imnln::printd("Start gaussian resampling inv");
        // Pt
        std::vector<complex_type> ptv = gaussian_resampling_Pt(tv, s, t);
        for (auto &i : ptv) {
            std::cout << i << std::endl;
        }
        std::cout << std::endl;

        //C
        std::vector<complex_type> csv = gaussian_resampling_C(ptv, s, t);
        for (auto &i : csv) {
            std::cout << i << std::endl;
        }
        std::cout << std::endl;

        // I'
        auto isv = gaussian_resampling_I(csv, s, t, a, p);
        for (auto &i : isv) {
            std::cout << i << std::endl;
        }
        std::cout << std::endl;

        // D'
        auto dsv = gaussian_resampling_D(isv, s, t, a);
        for (auto &i : dsv) {
            std::cout << i << std::endl;
        }
        std::cout << std::endl;

        // Ps-1
        auto pssv = gaussian_resampling_Psinv(dsv, s, t);
        auto complex_cmp = [](const complex_type &t1, const complex_type &t2) {
            return std::real(t1) < std::real(t2);
        };
        // std::sort(pssv.begin(), pssv.end(), complex_cmp);
        for (size_t i = 0; i < s; i++) {
            std::cout << std::round(std::real((pssv[i]*(complex_type)std::pow(2, 2*a*a)))) << "   ";
            std::cout << sv[i] << "    ";
            std::cout << (pssv[i]*(complex_type)std::pow(2, 2*a*a))/sv[i] << std::endl;
        }
        std::cout << std::endl;

        imnln::printd("End gaussian resampling inv");
        return 0;
    }

    /**
     * Given two intergers (as string) multiply them
     * using Schönhage-Strassen
     */
    int HVDH(std::string fname1, std::string fname2, std::string fout) {

        // Get first polynomial
        HVDHPolynomial<poly_type> p1;
        {
            std::vector<poly_type> v1;
            {
                std::string istr1;
                istr1 = imnln::read_integer(imnln::INTEGER_FILE_1);
                v1 = imnln::split<poly_type>(istr1);
                printd("Degree: " + std::to_string(v1.size()) + " S: " 
                        + std::to_string(imnln::params.S));
                v1.resize(imnln::params.S);
                agarwal_cooley(v1);
                istr1.erase(istr1.begin(), istr1.end());
            }
            p1.vec_to_poly(v1);
        }
        gaussian_resampling(p1);
        return 0;

        // Get second polynomial
        HVDHPolynomial<poly_type> p2;
        {
            std::vector<poly_type> v2;
            {
                std::string istr2;
                istr2 = imnln::read_integer(imnln::INTEGER_FILE_2);
                v2 = imnln::split<poly_type>(istr2);
                istr2.erase(istr2.begin(), istr2.end());
            }
            p2.vec_to_poly(v2);
        }

        imnln::print_ram_info();

        auto pr = p1.polymul(p2);

        // trim 0 coeffients at beginning
        size_t pr_len = pr.size()-1;
        while(pr_len > 0 && std::abs(std::round(pr[pr_len])) < imnln::EPSILON) {
            pr_len--;
        }
        // perform carrying
        uint64_t mod = std::pow(10, imnln::CHUCK_SIZE);
        for (size_t i = pr_len; i > 0 && i < ULLONG_MAX; i--) {
            uint64_t v = (uint64_t)std::round(pr[i]);
            pr[i-1] += (poly_type)(v / mod); // carry
            pr[i] = (poly_type)(v % mod);    // rest
        }

        // Build result string
        imnln::write_integer<poly_type>(fout, pr.get_coeffs(), pr_len);

        return 0;
    }

} // namespace imnln
