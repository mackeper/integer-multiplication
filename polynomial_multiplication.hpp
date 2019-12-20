// Authors: Marcus Östling, (& Tomas Möre, parts is from the course popup 2019)
#pragma once
#include <cmath>
#include <complex>
#include <valarray>
#include <algorithm>
#include <cstdio>
#include <vector>

namespace imnln {

    /**
     *  Class that stores the coeffients of a polynomial and using
     *  FFT can be multiplied with another Polynomial.
     */
    template <class T>
    class Polynomial {

        T* coeffs_;
        size_t size_;

        /**
         *  Internal Descrete Fourier Transform, solved by dividing and conquer.
         *  arr: valarray of complex coeffients.
         *  omegas: Complex point used for transformation. (Performance reason)
         */
        void DFT(std::valarray<std::complex<T>> &arr, std::valarray<std::complex<T>> &omegas) const {
            if(arr.size() <= 1) return;

            std::valarray<std::complex<T>> even = arr[std::slice(0, arr.size()/2, 2)];
            std::valarray<std::complex<T>> odd = arr[std::slice(1, arr.size()/2, 2)];
            std::valarray<std::complex<T>> new_omegas = omegas[std::slice(0, arr.size()/2, 2)];
            DFT(even, new_omegas);
            DFT(odd, new_omegas);
            for(size_t i = 0; i < arr.size()/2; i++) {
                auto tmp = omegas[i]*odd[i];
                arr[i] = even[i]+tmp;
                arr[i+arr.size()/2] = even[i]-tmp;
            }
        }

        public:
        explicit Polynomial(size_t size) {
            coeffs_ = new T[size];
            size_ = size;
        }

        explicit Polynomial(std::vector<T> v) {
            coeffs_ = new T[v.size()];
            size_ = v.size();
            for (size_t i = 0; i < v.size(); i++) 
                coeffs_[i] = v[i];
        }

        Polynomial(Polynomial<T> &&p) {
            this->size_ = p.size_;
            this->coeffs_ = p.coeffs_;
            p.coeffs_ = nullptr;
        }

        ~Polynomial() {
            delete[] coeffs_;
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

        /**
         *  Multiplies this polynomial with another, returns a new polynomial.
         */
        Polynomial<T> polymul(const Polynomial<T> &other) const {
            size_t size = (size_t)std::pow(2,
                    std::ceil(std::log2(std::max(size_, other.size()))) + 1);

            const std::complex<T> c = 0;
            std::valarray<std::complex<T>> p1(c, size);
            std::valarray<std::complex<T>> p2(c, size);

            // Init our complex valarrays
            for(size_t i = 0; i < size_; i++)
                p1[i] = std::complex<T>((*this)[i], 0);

            for(size_t i = 0; i < other.size(); i++)
                p2[i] = std::complex<T>(other[i], 0);

            // Init roots
            std::valarray<std::complex<T>> omegas(size);
            for(size_t i = 0; i < size; i++)
                omegas[i] = std::polar(1.0, -2.0 * M_PI * ((T)i/(T)size));

            // Run DFT
            DFT(p1, omegas);
            DFT(p2, omegas);

            // Dot multiplication
            for(size_t i = 0; i < size; i++)
                p1[i] *= p2[i];

            // Inverse DFT
            Polynomial p(size);
            for(size_t i = 0; i < size; i++)
                omegas[i] = std::polar(1.0, 2.0 * M_PI * ((T)i/(T)size));

            DFT(p1, omegas);

            for(size_t i = 0; i < size; i++) {
                p[i] = p1[i].real()/(T)size;
            }
            return p;
        }
    };

    /**
     * Given two intergers (as string) multiply them
     * using Schönhage-Strassen
     */
    std::string SSA(std::string istr1, std::string istr2) {
        Polynomial<double> p1(imnln::split<double>(istr1));
        Polynomial<double> p2(imnln::split<double>(istr2));
        auto pr = p1.polymul(p2);

        // trim 0 coeffients at beginning
        size_t pr_len = pr.size();
        while(pr_len > 0 && std::round(pr[pr_len]) == 0) {
            pr_len--;
        }

        // perform carrying
        uint64_t mod = std::pow(10, imnln::CHUCK_SIZE);
        for (size_t i = pr_len; i > 0 && i < ULLONG_MAX; i--) {
            uint64_t v = (uint64_t)std::round(pr[i]);
            pr[i-1] += v / mod; // carry
            pr[i] = v % mod;    // rest
        }

        // Build result string
        std::string result_str = "";
        for (size_t i = 0; i <= pr_len; i++) {
            uint64_t v = (uint64_t)std::round(pr[i]);
            std::string prs = std::to_string(v);
            if (i != 0)
                prs = std::string(imnln::CHUCK_SIZE - prs.length(), '0') + prs;
            result_str += prs;
        }

        return result_str;

    }

} // namespace imnln
