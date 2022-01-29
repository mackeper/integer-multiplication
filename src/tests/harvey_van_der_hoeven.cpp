#include <cstdint>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>

#include "util.hpp"
#include "types.hpp"
#include "catch.hpp"
#include "harvey_van_der_hoeven.hpp"

#define EPS 1e-15

TEST_CASE("to_unit_disc") {
    // Arrange
    std::vector<imnln::poly_type> vec;
    vec.push_back(0.0);
    vec.push_back(M_PI);
    vec.push_back(M_PI/2);
    vec.push_back(-M_PI/2);

    // Act
    auto vec2 = imnln::to_unit_disc(vec);

    // Assert
    REQUIRE(vec2[0].real() == Approx(1.0).margin(EPS));
    REQUIRE(vec2[0].imag() == Approx(0.0).margin(EPS));

    REQUIRE(vec2[1].real() == Approx(-1.0).margin(EPS));
    REQUIRE(vec2[1].imag() == Approx(0.0).margin(EPS));
    
    REQUIRE(vec2[2].real() == Approx(0.0).margin(EPS));
    REQUIRE(vec2[2].imag() == Approx(1.0).margin(EPS));

    REQUIRE(vec2[3].real() == Approx(0.0).margin(EPS));
    REQUIRE(vec2[3].imag() == Approx(-1.0).margin(EPS));
}