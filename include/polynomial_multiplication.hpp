// Authors: Marcus Östling, (& Tomas Möre, parts is from the course popup 2019)
#pragma once
#include <cmath>
#include <complex>
#include <valarray>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <fstream>

namespace imnln {

    typedef long double poly_type;

    /**
     *  Class that stores the coeffients of a polynomial and using
     *  FFT can be multiplied with another Polynomial.
     */
    template <class T>
    class Polynomial {

        T* coeffs_;
        size_t size_;
        bool deleted = true;

        /**
         *  Internal Descrete Fourier Transform, solved by dividing and conquer.
         *  arr: valarray of complex coeffients.
         *  omegas: Complex point used for transformation. (Performance reason)
         *
         *  Can probabilty be optimizied for memory
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
        Polynomial() {
            coeffs_ = new T[0];
            size_ = 0;
            deleted = false;
        }

        explicit Polynomial(size_t size) {
            coeffs_ = new T[size];
            size_ = size;
            deleted = false;
        }

        explicit Polynomial(std::vector<T> &v) {
            coeffs_ = new T[v.size()];
            size_ = v.size();
            for (size_t i = 0; i < v.size(); i++) {
                coeffs_[i] = v[i];
            }
            deleted = false;
        }

        Polynomial(Polynomial<T> &&p) {
            this->size_ = p.size_;
            this->coeffs_ = p.coeffs_;
            p.coeffs_ = nullptr;
            deleted = false;
        }

        ~Polynomial() {
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
        Polynomial<T> polymul(Polynomial<T> &other) {
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
            Polynomial p(size);
            imnln::print_ram_info();
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
    int SSA(std::string fname1, std::string fname2, std::string fout) {

        // Get first polynomial
        Polynomial<poly_type> p1;
        {
            std::vector<poly_type> v1;
            {
                std::string istr1;
                istr1 = imnln::read_integer(imnln::INTEGER_FILE_1);
                v1 = imnln::split(istr1, imnln::CHUCK_SIZE);
                istr1.erase(istr1.begin(), istr1.end());
            }
            p1.vec_to_poly(v1);
        }

        // Get second polynomial
        Polynomial<poly_type> p2;
        {
            std::vector<poly_type> v2;
            {
                std::string istr2;
                istr2 = imnln::read_integer(imnln::INTEGER_FILE_2);
                v2 = imnln::split(istr2, imnln::CHUCK_SIZE);
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
        imnln::write_integer(fout, pr.get_coeffs(), pr_len, imnln::CHUCK_SIZE);

        return 0;
    }

} // namespace imnln
