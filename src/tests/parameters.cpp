#include <iostream>
#include <stdexcept>
#include "types.hpp"
#include "catch.hpp"
#include "parameters.hpp"

#define EPS 1e-15

TEST_CASE("unique_power_of_two") {
    REQUIRE(imnln::unique_power_of_two(3, 5) == 4);
    REQUIRE(imnln::unique_power_of_two(50, 80) == 64);
    REQUIRE(imnln::unique_power_of_two(0, 1) == 0);
}

TEST_CASE("is_prime") {
    REQUIRE(imnln::is_prime(10019689));
    REQUIRE_FALSE(imnln::is_prime(10049681));
    REQUIRE_FALSE(imnln::is_prime(0));
    REQUIRE_FALSE(imnln::is_prime(1));
    REQUIRE(imnln::is_prime(2));
    REQUIRE(imnln::is_prime(3));
    REQUIRE_FALSE(imnln::is_prime(4));
}

TEST_CASE("find_prime_under") {
    REQUIRE(imnln::find_prime_under(10, 0) == 7);
    REQUIRE(imnln::find_prime_under(10, 1) == 5);
    REQUIRE(imnln::find_prime_under(10049681, 0) == 10049671);
    REQUIRE(imnln::find_prime_under(10049681, 1) == 10049657);
    REQUIRE(imnln::find_prime_under(10, 100) == 0);
}