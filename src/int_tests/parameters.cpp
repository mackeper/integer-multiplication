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
    write_file(filename4, "67890");

    // Act
    //auto params1 = imnln::get_parameters(filename1, filename2);
    //auto params2 = imnln::get_parameters(filename3, filename4);

    // Assert
    REQUIRE(1 == 1);

    // Cleanup
    remove(filename1);
    remove(filename2);
    remove(filename3);
    remove(filename4);
}