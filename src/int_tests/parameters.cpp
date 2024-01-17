#include <iostream>
#include <fstream>
#include <stdexcept>
#include "types.hpp"
#include "catch.hpp"
#include "parameters.hpp"

#define EPS 1e-15

TEST_CASE("get_parameters") {
    // Arrange
    auto const filename1 = "test_get_parameters1.txt";
    auto const filename2 = "test_get_parameters2.txt";
    auto const filename3 = "test_get_parameters3.txt";
    auto const filename4 = "test_get_parameters4.txt";

    auto write_file = [] (std::string filename, std::string content) {
        std::ofstream f;
        f.open(filename);
        f << content;
        f.close();
    };

    write_file(filename1, "1");
    write_file(filename2, "1");
    write_file(filename3, "12345");
    write_file(filename4, "6789");

    // Act
    auto params1 = imnln::get_parameters(filename1, filename2);
    auto params2 = imnln::get_parameters(filename3, filename4);

    // Cleanup
    remove(filename1);
    remove(filename2);
    remove(filename3);
    remove(filename4);

    // Assert
    REQUIRE(params1.d == 2);
    REQUIRE(params1.n == 4);
    REQUIRE(params1.b == 2);
    REQUIRE(params1.p == 12);
    REQUIRE(params1.alpha == 4);
    REQUIRE(params1.gamma == 64);
    REQUIRE(params1.T == 8);
    REQUIRE(params1.r == 4);
    REQUIRE(params1.t[0] == 2);
    REQUIRE(params1.t[1] == 4);
    REQUIRE(params1.S == 4);

    REQUIRE(params2.d == 2);
    REQUIRE(params2.n == 20);
    REQUIRE(params2.b == 5);
    REQUIRE(params2.p == 30);
    REQUIRE(params2.alpha == 4);
    REQUIRE(params2.gamma == 64);
    REQUIRE(params2.T == 16);
    REQUIRE(params2.r == 4);
    REQUIRE(params2.t[0] == 4);
    REQUIRE(params2.t[1] == 4);
}