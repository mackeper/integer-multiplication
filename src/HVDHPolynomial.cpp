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
#include <iostream>
#include <stdexcept>
#include "types.hpp"
#include "agarwal_cooley.hpp"
#include "HVDHPolynomial.hpp"

imnln::HVDHPolynomial::HVDHPolynomial() {
    size_ = 0;
    coeffs_ = new poly_type[size_];
    deleted = false;
}

imnln::HVDHPolynomial::HVDHPolynomial(size_t size) {
    coeffs_ = new poly_type[size];
    size_ = size;
    deleted = false;
}

imnln::HVDHPolynomial::HVDHPolynomial(std::vector<poly_type> &v) {
    coeffs_ = new poly_type[v.size()];
    size_ = v.size();
    for (size_t i = 0; i < v.size(); i++) {
        coeffs_[i] = v[i];
    }
    deleted = false;
}

imnln::HVDHPolynomial::HVDHPolynomial(const HVDHPolynomial &p) {
    this->size_ = p.size_;
    this->coeffs_ = new poly_type[p.size_];
    for (size_t i = 0; i < p.size_; i++) {
        coeffs_[i] = p.coeffs_[i];
    }
    deleted = false;
}

imnln::HVDHPolynomial::HVDHPolynomial(HVDHPolynomial &&p) {
    this->size_ = p.size_;
    this->coeffs_ = p.coeffs_;
    p.coeffs_ = nullptr;
    p.size_ = 0;
    deleted = false;
}

imnln::HVDHPolynomial::~HVDHPolynomial() {
    if (!deleted) {
        delete[] coeffs_;
        size_ = 0;
        deleted = true;
    }
}

void imnln::HVDHPolynomial::erase() {
    if (!deleted) {
        delete[] coeffs_;
        size_ = 0;
        deleted = true;
    }
}

void imnln::HVDHPolynomial::vec_to_poly(std::vector<poly_type> &v) {
    this->erase();
    coeffs_ = new poly_type[v.size()];
    size_ = v.size();
    for (size_t i = 0; i < v.size(); i++) {
        coeffs_[i] = v[i];
    }
    deleted = false;
}

size_t imnln::HVDHPolynomial::size() const {
    return size_;
}

imnln::poly_type& imnln::HVDHPolynomial::operator[](size_t i) {
    if (i >= size_)
        throw std::out_of_range("Index out of range.");
    return coeffs_[i];
}

imnln::poly_type imnln::HVDHPolynomial::operator[](size_t i) const {
    if (i >= size_)
        throw std::out_of_range("Index out of range.");
    return coeffs_[i];
}

imnln::poly_type* imnln::HVDHPolynomial::get_coeffs() {
    return coeffs_;
}