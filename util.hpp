#ifndef INTEGER_MULTIPLICATION_UTIL
#define INTEGER_MULTIPLICATION_UTIL
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <climits>

#include "parameters.hpp"

// process memory
#include "sys/types.h"
#include "sys/sysinfo.h"
#include <iomanip>

#define DEBUG 0

namespace imnln {
    auto start = std::chrono::high_resolution_clock::now();
    auto stop = std::chrono::high_resolution_clock::now();

    void timer_start() {
        start = std::chrono::high_resolution_clock::now();
    }

    void timer_stop(std::string text = "") {
        stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "\t" << text << duration.count() << std::endl;

    }

    void printd(std::string s) {
        if (DEBUG)
            std::cout << s <<std::endl;
    }

    void generate_integer(std::string fname, size_t length) {
        printd("generate_integer " + fname + " start");
        std::ofstream f;
        f.open(fname);
        for (size_t i = 0; i < length; i++) {
            f << std::to_string(rand()%10);
        }
        f.close();
        printd("generate_integer " + fname + " end");
    }

    std::string read_integer(std::string fname) {
        printd("read_integer " + fname + " start");
        std::string in_int;
        std::ifstream f(fname);
        std::string buf;
        if(f.is_open()) {
            while(f >> buf) {
                in_int.append(buf);
            }
        }
        f.close();
        printd("read_integer " + fname + " end");
        return in_int;
    }

    bool compare_integer(std::string fname1, std::string fname2) {
        std::string tmp1 = read_integer(fname1);
        std::string tmp2 = read_integer(fname2);
        return tmp1.compare(tmp2) == 0;
    }

    template <class T>
    int write_integer(std::string fname,  T* value, size_t len) {
        printd("write_integer " + fname + " start");

        std::ofstream f;
        f.open(fname);
        for (size_t i = 0; i <= len; i++) {
            uint64_t v = (uint64_t)std::round(value[i]);
            std::string prs = std::to_string(v);
            if (i != 0)
                prs = std::string(imnln::CHUCK_SIZE - prs.length(), '0') + prs;
            f << prs;
        }
        f.close();

        printd("write_integer " + fname + " end");
        return 0;
    }

    int write_integer(std::string fname,  char * value) {
        printd("write_integer " + fname + " start");

        std::ofstream f;
        f.open(fname);
        f << value;
        f.close();

        printd("write_integer " + fname + " end");
        return 0;
    }


    template <class T>
    T str_to_num(std::string &istr) {
        T i = 0;

        for (char c : istr) {
            i = i*10 + (int)c - 48;
        }

        return i;
    }

    template <class T>
    std::vector<T> split(std::string &istr) {
        printd("split start");
        std::vector<T> v;
        size_t i = 0;
        size_t j = (imnln::CHUCK_SIZE - istr.size() % imnln::CHUCK_SIZE) % imnln::CHUCK_SIZE;
        std::string tmp_str = "";
        while(i < istr.size()) {
            if (j == imnln::CHUCK_SIZE) {
                v.push_back(str_to_num<T>(tmp_str)); 
                tmp_str = "";
                j = 0;
            }
            tmp_str.push_back(istr[i]);
            i++;
            j++;
        }
        if (tmp_str.compare("") != 0)
                v.push_back(str_to_num<T>(tmp_str)); 

        // std::cout << std::fixed;
        // for (auto e : v)
        //     std::cout << e << std::endl;
        // std::cout << std::scientific;
        printd("split end");
        return v;
    }

    template <class T>
    std::string concat(std::vector<T> &v) {
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
    void print_ram_info(){
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

}

#endif // INTEGER_MULTIPLICATION_UTIL

