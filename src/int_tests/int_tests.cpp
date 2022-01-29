#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#define EPS 1e-15

TEST_CASE("main") {
    // Arrange
    auto a = 1;
    auto b = 2;

    // Act
    auto c = a + b;

    // Assert
    REQUIRE(c == Approx(3.0).margin(EPS));
}