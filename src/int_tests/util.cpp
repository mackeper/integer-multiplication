#include <iostream>
#include <fstream>
#include <sstream>
#include "types.hpp"
#include "catch.hpp"
#include "util.hpp"

#define EPS 1e-15

TEST_CASE("generate_integer") {
    // Arrange
    auto const filename1 = "tmp_test_integer1.txt";
    auto const filename2 = "tmp_test_integer2.txt";
    auto const filename3SameAs1 = "tmp_test_integer3.txt";
    auto const filename4DifferentFrom1 = "tmp_test_integer4.txt";

    // Act
    imnln::generate_integer(filename1, 20);
    imnln::generate_integer(filename2, 35, 3);
    imnln::generate_integer(filename3SameAs1, 20);
    imnln::generate_integer(filename4DifferentFrom1, 20, 2);

    // Assert
    auto read_file = [](std::string filename) {
        std::string in_int;
        std::ifstream f(filename);
        std::string buf;
        if(f.is_open()) {
            while(f >> buf) {
                in_int.append(buf);
            }
        }
        f.close();
        return in_int;
    };

    REQUIRE(read_file(filename1) == "36753562912709360626");
    REQUIRE(read_file(filename2) == "65805026148439214515528837578844646");
    REQUIRE(read_file(filename3SameAs1) == "36753562912709360626");
    REQUIRE(read_file(filename4DifferentFrom1) == "09851847572993794370");
    REQUIRE(read_file(filename1) == read_file(filename3SameAs1));
    REQUIRE(read_file(filename1) != read_file(filename4DifferentFrom1));

    // Cleanup
    remove(filename1);
    remove(filename2);
    remove(filename3SameAs1);
    remove(filename4DifferentFrom1);
}

TEST_CASE("read_integer") {
    // Arrange
    auto const filename = "test_read_integer.txt";
    std::ofstream f;
    f.open(filename);
    for (size_t i = 0; i < 11; i++) {
        f << std::to_string(i);
    }
    f.close();

    // Act
    auto in_int = imnln::read_integer(filename);

    // Assert
    REQUIRE(in_int == "012345678910");

    // Cleanup
    remove(filename);
}

TEST_CASE("compare_integer") {
    // Arrange
    auto const filename1 = "test_compare_integer1.txt";
    auto const filename2 = "test_compare_integer2.txt";
    auto const filename3 = "test_compare_integer3.txt";
    auto const filename4 = "test_compare_integer4.txt";

    auto write_file = [] (std::string filename, std::string content) {
        std::ofstream f;
        f.open(filename);
        f << content;
        f.close();
    };

    write_file(filename1, "012345678910");
    write_file(filename2, "012345678910");
    write_file(filename3, "012345678911");
    write_file(filename4, "0123456789111");

    // Act
    auto result1 = imnln::compare_integer(filename1, filename2);
    auto result2 = imnln::compare_integer(filename1, filename3);
    auto result3 = imnln::compare_integer(filename2, filename3);
    auto result4 = imnln::compare_integer(filename3, filename4);

    // Assert
    REQUIRE(result1);
    REQUIRE(!result2);
    REQUIRE(!result3);
    REQUIRE(!result4);

    // Cleanup
    remove(filename1);
    remove(filename2);
    remove(filename3);
    remove(filename4);
}

TEST_CASE("write_integer poly_tyme") {
    // Arrange
    auto const size1 = 9;
    auto const chuck_size1 = 4;
    auto const filename1 = "test_write_integer1.txt";
    imnln::poly_type* poly_type1 = new imnln::poly_type[size1]
        {1, 2, 3, 4, 5, 6, 7, 8, 9};

    auto const size2 = 16;
    auto const chuck_size2 = 10;
    auto const filename2 = "test_write_integer2.txt";
    imnln::poly_type* poly_type2 = new imnln::poly_type[size2]
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

    // Act
    imnln::write_integer(filename1, poly_type1, size1, chuck_size1);
    imnln::write_integer(filename2, poly_type2, size2, chuck_size2);

    // Assert
    auto read_file = [](std::string filename) {
        std::string in_int;
        std::ifstream f(filename);
        std::string buf;
        if(f.is_open()) {
            while(f >> buf) {
                in_int.append(buf);
            }
        }
        f.close();
        return in_int;
    };

    auto in_int1 = read_file(filename1);
    auto result1 = false;
    for (size_t i = 0; i < size1; i++)
        result1 |= (int)(in_int1[i*chuck_size1]-'0') != (int)poly_type1[i];
    REQUIRE_FALSE(result1);

    auto in_int2 = read_file(filename2);
    auto result2 = false;
    for (size_t i = 0; i < size2; i++)
        result2 |= (int)(in_int2[i*chuck_size2]-'0') != (int)poly_type2[i]%10;
    REQUIRE_FALSE(result2);

    // Cleanup
    remove(filename1);
    remove(filename2);
}

TEST_CASE("write_integer char*") {
    // Arrange
    auto const filename1 = "test_write_integer_char1.txt";
    auto const filename2 = "test_write_integer_char2.txt";
    char* int1 = new char[10] {'1', '2', '3', '4', '5', '6', '7', '8', '9', 0};
    char* int2 = new char[12]
        {'1', '2', '3', '4', '5', '6', '7', '8', '9', '1', '0', 0};

    // Act
    imnln::write_integer(filename1, int1);
    imnln::write_integer(filename2, int2);

    // Assert
    auto read_file = [](std::string filename) {
        std::string in_int;
        std::ifstream f(filename);
        std::string buf;
        if(f.is_open()) {
            while(f >> buf) {
                in_int.append(buf);
            }
        }
        f.close();
        return in_int;
    };

    auto in_int1 = read_file(filename1);
    auto in_int2 = read_file(filename2);
    REQUIRE(in_int1 == "123456789");
    REQUIRE(in_int2 == "12345678910");

    // Cleanup
    remove(filename1);
    remove(filename2);
}