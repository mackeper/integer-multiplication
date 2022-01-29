#include <iostream>
#include <fstream>
#include <sstream>
#include "types.hpp"
#include "catch.hpp"
#include "util.hpp"

#define EPS 1e-15

TEST_CASE("str_to_num") {
    // Arrange
    std::string str1 = "123456789";
    const std::string str2 = "123456789101112";

    // Act
    auto num1 = imnln::str_to_num(str1);
    auto const num2 = imnln::str_to_num(str2);

    // Assert
    REQUIRE(num1 == imnln::poly_type(123456789));
    REQUIRE(num2 == imnln::poly_type(123456789101112));
}

TEST_CASE("split") {
    // Arrange
    std::string str1 = "123456789";
    const std::string str2 = "12345678910111213141516";

    // Act
    auto num1 = imnln::split(str1, 4);
    auto const num2 = imnln::split(str2, 4);

    // Assert
    const std::vector<imnln::poly_type> expected1
        {1, 2345, 6789};
    auto result1 = false;
    for (size_t i = 0; i < num1.size(); i++)
        result1 |= expected1[i] != num1[i];
    REQUIRE_FALSE(result1);

    const std::vector<imnln::poly_type> expected2
        {123, 4567, 8910, 1112, 1314, 1516};
    auto result2 = false;
    for (size_t i = 0; i < num2.size(); i++)
        result2 |= expected2[i] != num2[i];
    REQUIRE_FALSE(result2);
}

TEST_CASE("concat") {
    // Arrange
    std::vector<imnln::poly_type> vec1 {1, 2, 3};
    std::vector<imnln::poly_type> vec2;
    std::vector<imnln::poly_type> vec3 {1};

    // Act
    auto result1 = imnln::concat(vec1);
    auto result2 = imnln::concat(vec2);
    auto result3 = imnln::concat(vec3);

    // Assert
    REQUIRE(result1 == "123");
    REQUIRE(result2 == "");
    REQUIRE(result3 == "1");
}

TEST_CASE("transpose") {
    // Arrange
    std::vector<imnln::poly_type> vec1 {1};
    std::vector<imnln::poly_type> vec2 {1, 2};
    std::vector<imnln::poly_type> vec3 {1, 2, 3, 4};

    // Act
    imnln::transpose(vec1, 1, 1);
    imnln::transpose(vec2, 1, 2);
    imnln::transpose(vec3, 2, 2);

    // Assert
    REQUIRE(vec1[0] == 1);

    REQUIRE(vec2[0] == 1);
    REQUIRE(vec2[1] == 2);

    REQUIRE(vec3[0] == 1);
    REQUIRE(vec3[1] == 3);
    REQUIRE(vec3[2] == 2);
    REQUIRE(vec3[3] == 4);
}

TEST_CASE("mod_inverse") {
    REQUIRE(imnln::mod_inverse(123, 4567) == 854);
    REQUIRE(imnln::mod_inverse(2, 3) == 2);
    REQUIRE(imnln::mod_inverse(2, 2) == 0);
}