#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <climits>
#include <bitset>
#include <chrono>
#include <iostream>
#include <fstream>

// process memory
#include "sys/types.h"
#include "sys/sysinfo.h"
#include <iomanip>

#include "parameters.hpp"
#include "util.hpp"

auto static start = std::chrono::high_resolution_clock::now();
auto static stop = std::chrono::high_resolution_clock::now();

void imnln::timer_start() {
    start = std::chrono::high_resolution_clock::now();
}

void imnln::timer_stop(const std::string text) {
    stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "\ttimer: " << text << duration.count() << std::endl;

}

void imnln::printd(const std::string s) {
    #ifndef DEBUG
    #define DEBUG false
    #endif
    if (DEBUG)
        std::cout << s <<std::endl;
}

void imnln::generate_integer(const std::string fname, const size_t length, const int seed /*= 0*/) {
    srand(seed);
    std::ofstream f;
    f.open(fname);
    for (size_t i = 0; i < length; i++) {
        f << std::to_string(rand()%10);
    }
    f.close();
}

std::string imnln::read_integer(const std::string fname) {
    //printd("read_integer " + fname + " start");
    std::string in_int;
    std::ifstream f(fname);
    std::string buf;
    if(f.is_open()) {
        while(f >> buf) {
            in_int.append(buf);
        }
    }
    f.close();
    //printd("read_integer " + fname + " end");
    return in_int;
}

bool imnln::compare_integer(const std::string fname1, const std::string fname2) {
    std::string tmp1 = read_integer(fname1);
    std::string tmp2 = read_integer(fname2);
    return tmp1.compare(tmp2) == 0;
}

int imnln::write_integer(const std::string fname, const poly_type* value, const size_t len, const size_t chuck_size) {
    //printd("write_integer " + fname + " start");

    std::ofstream f;
    f.open(fname);
    for (size_t i = 0; i <= len; i++) {
        uint64_t v = (uint64_t)std::round(value[i]);
        std::string prs = std::to_string(v);
        if (i != 0)
            prs = std::string(chuck_size - prs.length(), '0') + prs;
        f << prs;
    }
    f.close();

    //printd("write_integer " + fname + " end");
    return 0;
}

int imnln::write_integer(const std::string fname, const char * value) {
    //printd("write_integer " + fname + " start");

    std::ofstream f;
    f.open(fname);
    f << value;
    f.close();

    //printd("write_integer " + fname + " end");
    return 0;
}


imnln::poly_type imnln::str_to_num(const std::string &istr) {
    imnln::poly_type i = 0;

    for (char c : istr) {
        i = i*10 + (int)c - 48;
    }

    return i;
}

std::vector<imnln::poly_type> imnln::split(const std::string &istr, const size_t chuck_size) {
    //printd("split start");
    std::vector<imnln::poly_type> v;
    size_t i = 0;
    size_t j = (chuck_size - istr.size() % chuck_size) % chuck_size;
    std::string tmp_str = "";
    while(i < istr.size()) {
        if (j == chuck_size) {
            v.push_back(str_to_num(tmp_str)); 
            tmp_str = "";
            j = 0;
        }
        tmp_str.push_back(istr[i]);
        i++;
        j++;
    }
    if (tmp_str.compare("") != 0)
            v.push_back(str_to_num(tmp_str)); 

    // std::cout << std::fixed;
    // for (auto e : v)
    //     std::cout << e << std::endl;
    // std::cout << std::scientific;
    //printd("split end");
    return v;
}

std::string imnln::concat(const std::vector<poly_type> &v) {
    std::string is = "";
    for(auto it = v.begin(); it != v.end(); it++) {
        std::string cs = "";
        cs += std::to_string((uint64_t)std::round(*it));
        // while (cs.size() < imnln::CHUCK_SIZE) {
        //     cs.insert(cs.begin(), '0');
        // }
        is += cs;
    }
    return is;
}

// source: https://gist.github.com/thirdwing/da4621eb163a886a03c5
void imnln::print_ram_info(){
    struct sysinfo memInfo;
    sysinfo (&memInfo);
    // double totalVirtualMem = memInfo.totalram;
    double virtualMemUsed = (double)memInfo.totalram - (double)memInfo.freeram;

    virtualMemUsed /= 1024;
    virtualMemUsed /= 1024;
    virtualMemUsed /= 1024;
    
    if (DEBUG) {
        std::cout << std::fixed;
        std::cout << std::setprecision(3);
        std::cout << "Used ram: " << virtualMemUsed << " G" << std::endl;
    }

}

// Non-square matrix transpose of matrix of size r x c and base address A 
// https://www.geeksforgeeks.org/inplace-m-x-n-size-matrix-transpose/
void imnln::transpose(std::vector<imnln::poly_type> &v, const int r, const int c) { 
    int size = r*c - 1; 
    poly_type t; // holds element to be replaced, eventually becomes next element to move 
    int next; // location of 't' to be moved 
    int cycleBegin; // holds start of cycle 
    int i; // iterator 
    // hash to mark moved elements 
    std::vector<bool> b(v.size());
    b[0] = b[size] = 1; 
    i = 1; // Note that A[0] and A[size-1] won't move 
    while (i < size) { 
        cycleBegin = i; 
        t = v[i]; 
        do { 
            // Input matrix [r x c] 
            // Output matrix  
            // i_new = (i*r)%(N-1) 
            next = (i*r)%size; 
            std::swap(v[next], t); 
            b[i] = 1; 
            i = next; 
        } 
        while (i != cycleBegin); 
    
        // Get Next Move (what about querying random location?) 
        for (i = 1; i < size && b[i]; i++) 
            ; 
    } 
}

// https://www.geeksforgeeks.org/multiplicative-inverse-under-modulo-m/
// Inverse of a modulo m
uint64_t imnln::mod_inverse(const uint64_t a, const uint64_t m) { 
    u_int64_t a2 = a;
    a2 = a2%m; 
    for (uint64_t x = 1; x<m; x++) { 
        if ((a2*x) % m == 1) {
            return x; 
        }
    }
    return 0;
}