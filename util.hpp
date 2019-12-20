#ifndef INTEGER_MULTIPLICATION_UTIL
#define INTEGER_MULTIPLICATION_UTIL
#include <cstdint>
#include <string>
#include <vector>
#include <cmath>
#include <climits>

#include "parameters.hpp"

namespace imnln {
    template <class T>
    T str_to_num(std::string istr) {
        T i = 0;

        for (char c : istr) {
            i = i*10 + (int)c - 48;
        }

        return i;
    }

    template <class T>
    std::vector<T> split(std::string istr) {
        std::vector<T> v;
        size_t i = istr.size()-1;
        size_t last_i = istr.size();
        while(i >= 0 && i < ULLONG_MAX) {
            // TODO insert slow?
            std::string tmp_str = "";
            for (int j = 0; j < imnln::CHUCK_SIZE && i >= 0 && i < ULLONG_MAX; j++, i--) {
                tmp_str.insert(tmp_str.begin(), istr[i]);
            }
            v.insert(v.begin(), str_to_num<T>(tmp_str)); 
        }

        return v;
    }

    template <class T>
    std::string concat(std::vector<T> v) {
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

}

#endif // INTEGER_MULTIPLICATION_UTIL

