#include <cstdint>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>

#include "util.hpp"
#include "types.hpp"
#include "catch.hpp"
#include "gaussian_resampling.hpp"

#define EPS 1e-15

TEST_CASE("gaussian_resampling_S") {
    // Arrange
    std::vector<imnln::complex_type> vec(10, 1);

    // Act
    auto vec2 = imnln::gaussian_resampling_S(vec, 10, 13, 2, 100);

    // Assert
    //REQUIRE(vec2[0].real() == Approx(1.0).margin(EPS));
    //REQUIRE(vec2[0].imag() == Approx(0.0).margin(EPS));

    //REQUIRE(vec2[1].real() == Approx(1.0).margin(EPS));
    //REQUIRE(vec2[1].imag() == Approx(0.0).margin(EPS));
    
    //REQUIRE(vec2[2].real() == Approx(1.0).margin(EPS));
    //REQUIRE(vec2[2].imag() == Approx(0.0).margin(EPS));

    //REQUIRE(vec2[3].real() == Approx(1.0).margin(EPS));
    //REQUIRE(vec2[3].imag() == Approx(0.0).margin(EPS));
}