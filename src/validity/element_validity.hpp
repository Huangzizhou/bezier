#pragma once
#include <array>
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"

namespace element_validity {
	Validity isValid(uint n, uint s, uint p, bool t, const std::vector<fp_t> &cp);
}
