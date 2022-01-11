#include <iostream>
#include <fstream>
#include <set>
#include <tuple>
#include <algorithm>
#include "util.hpp"
#include "types.hpp"
#include "parameters.hpp"
#include "agarwal_cooley.hpp"

int imnln::agarwal_cooley(std::vector<poly_type> &v, parameters params) {
	std::vector<uint64_t> s = params.s;
	std::vector<uint64_t> t = params.t;
	uint64_t S = params.S;
	uint64_t T = params.T;

	std::vector<uint64_t> v2 = {0, 1, 2, 3, 4, 4, 6, 7, 8, 9};
	std::vector<uint64_t> invs(s.size());
	for (uint64_t i = 0; i < s.size(); i++) {
		invs[i] = imnln::modInverse(t[i], s[i]);
		std::cout << "invs: " << invs[i] << std::endl;
	}

	std::vector<std::tuple<uint64_t, uint64_t>> v3;
	// size_t N = 10;
	// size_t r1 = 2;
	// size_t r2 = 5;
	// size_t c = 2;


	// uint64_t s1 = modInverse(r2, r1) * r2;
	// uint64_t s2 = modInverse(r1, r2) * r1;
	// std::cout << "s1: " << s1 << " s2: " << s2 << "\n";
	// for (auto i : v2) {
	//     uint64_t i1 = i % r1;
	//     uint64_t i2 = i % r2;
	//     for (uint64_t x = 0; x < s.size(); x++) {
	//
	//     }
	//     uint64_t inv = (i1*s1 + i2*s2) % N;
	//     for (auto i : s) {
	//
	//     }
	//
	//     std::cout << "i: " << i << " i1: " << i1 << " i2: " << i2 <<
	//         " inverse: " << inv << "\n";
	//     v3.push_back(std::make_tuple(i1, i2));
	//     assert(i == inv);
	// }
	//
	// for (uint64_t k2 = 0; k2 < r2; k2++) {
	//     for (uint64_t k1 = 0; k1 < r1; k1++) {
	//
	//     }
	// }

	// return 0;
	//std::vector<int> v2 = {1,2,3,4,5,6};
	// for (size_t i = 0; i < n; i++) {
	//     for (size_t j = 0; j < m; j++) {
	//         std::cout << v2[i*m+j] << " ";
	//     }
	//     std::cout << std::endl;
	// }
	// std::cout << std::endl;
	// imnln::transpose(v2, n, m);
	// for (size_t i = 0; i < m; i++) {
	//     for (size_t j = 0; j < n; j++) {
	//         std::cout << v2[i*n+j] << " ";
	//     }
	//     std::cout << std::endl;
	// }
	
	// std::cout << "Y: " << std::endl;
	// for (size_t i = 0; i < n; i++)  {
	//     for (size_t j = 0; j < m; j++)  {
	//         // std::cout << v2[(((i-j*c)%m)*m+j)*c] << std::endl;
	//         std::cout << v2[((i-j*c)%m)*m+j] << std::endl;
	//     }
	// }
	//
	// return 0;
	//
	std::vector<std::tuple<uint64_t, poly_type>> label_and_sort(S);
	std::set<uint64_t> st;
	for (size_t i = 0; i < S; i++) {
		label_and_sort[i] = std::tuple<poly_type, poly_type>((T*i) % S, v[i]);
		st.insert((T*i)%S);
	}

	std::sort(label_and_sort.begin(), label_and_sort.end());

	for (size_t i = 0; i < S; i++) {
		// std::cout << std::get<0>(label_and_sort[i]) << ", "
		//     << std::get<1>(label_and_sort[i]) << std::endl;
		if (invs[0] * i < S) {
		std:: cout << invs[0] * i << std::endl;
		// std::cout << std::get<1>(label_and_sort[invs[0]*i]) << std::endl;
		}
	}
	return 0;
}