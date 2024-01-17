#ifndef IMNLN_AGARWAL_COOLEY
#define IMNLN_AGARWAL_COOLEY
#include <vector>
#include "types.hpp"
#include "parameters.hpp"

namespace imnln {
	int agarwal_cooley(std::vector<poly_type> &v, parameters params);
}

#endif