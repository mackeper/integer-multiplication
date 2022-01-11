#ifndef INTEGER_MULTIPLICATION_HVDHPOLYNOMIAL
#define INTEGER_MULTIPLICATION_HVDHPOLYNOMIAL
#include <vector>
#include "types.hpp"
#include "agarwal_cooley.hpp"

namespace imnln {
    /**
     *  Class that stores the coeffients of a polynomial and using
     *  FFT can be multiplied with another HVDHPolynomial.
     */
    class HVDHPolynomial {
        poly_type* coeffs_;
        size_t size_;
        bool deleted = true;

        public:
        HVDHPolynomial();
        explicit HVDHPolynomial(size_t size);
        explicit HVDHPolynomial(std::vector<poly_type> &v);
        HVDHPolynomial(HVDHPolynomial &&p);
        ~HVDHPolynomial();

        void erase();
        void vec_to_poly(std::vector<poly_type> &v);
        size_t size() const;
        poly_type& operator[](size_t i);
        poly_type operator[](size_t i) const;
        poly_type* get_coeffs();
    };
}
#endif