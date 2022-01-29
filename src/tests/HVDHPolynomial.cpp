#include <iostream>
#include <stdexcept>
#include "types.hpp"
#include "catch.hpp"
#include "HVDHPolynomial.hpp"

#define EPS 1e-15

TEST_CASE("Constructor default") {
    // Arrange
    imnln::HVDHPolynomial p;

    // Act
    auto size = p.size();

    // Assert
    REQUIRE(size == 0);
}

TEST_CASE("Constructor size") {
    // Arrange
    imnln::HVDHPolynomial p(3);

    // Act
    auto size = p.size();

    // Assert
    REQUIRE(size == 3);
    REQUIRE_NOTHROW(p[2]);
    REQUIRE_THROWS(p[3]);
}

TEST_CASE("Constructor vector") {
    // Arrange
    std::vector<imnln::poly_type> vec {1, 2, 3};
    imnln::HVDHPolynomial p(vec);

    // Act
    auto size = p.size();

    // Assert
    REQUIRE(size == 3);
    REQUIRE(p[2] == 3);
    REQUIRE_THROWS(p[3]);
}

TEST_CASE("Constructor copy") {
    // Arrange
    std::vector<imnln::poly_type> vec {1, 2, 3};
    imnln::HVDHPolynomial p1(vec);
    imnln::HVDHPolynomial p2(p1);

    // Act
    auto size = p2.size();

    // Assert
    REQUIRE(size == 3);
    REQUIRE(p2[2] == 3);
    REQUIRE_THROWS(p2[3]);
}

TEST_CASE("Constructor move") {
    // Arrange
    std::vector<imnln::poly_type> vec {2, 3, 4};
    imnln::HVDHPolynomial p1(vec);
    imnln::HVDHPolynomial p2(std::move(p1));

    // Act
    auto size = p2.size();

    // Assert
    REQUIRE(size == 3);
    REQUIRE(p2[2] == 4);
    REQUIRE_THROWS(p1[2]);
    REQUIRE_THROWS(p2[3]);
}
TEST_CASE("erase") {
    // Arrange
    std::vector<imnln::poly_type> vec {2, 3, 4};
    imnln::HVDHPolynomial p(vec);

    // Act
    p.erase();

    // Assert
    REQUIRE(p.size() == 0);
}

TEST_CASE("vec_to_poly") {
    // Arrange
    std::vector<imnln::poly_type> vec {1, 2, 3};
    imnln::HVDHPolynomial p;
    assert(p.size() == 0);

    // Act
    p.vec_to_poly(vec);

    // Assert
    REQUIRE(p.size() == 3);
    REQUIRE(p[2] == 3);
    REQUIRE_THROWS(p[3]);
}

TEST_CASE("get_coeffs") {
    // Arrange
    std::vector<imnln::poly_type> vec {2, 3, 4};
    imnln::poly_type* coeffs;
    imnln::HVDHPolynomial p(vec);

    // Act
    coeffs = p.get_coeffs();

    // Assert
    REQUIRE(coeffs[0] == 2);
    REQUIRE(coeffs[1] == 3);
    REQUIRE(coeffs[2] == 4);
}
