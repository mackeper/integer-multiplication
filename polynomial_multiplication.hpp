// Authors: Marcus Östling, Tomas Möre 2019
#pragma once
#include <cmath>
#include <complex>
#include <valarray>
#include <algorithm>
#include <cstdio>

namespace popup {

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
            size_t size = std::pow(2, std::ceil(std::log2(std::max(size_, other.size()))) + 1);

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
} // namespace popup
