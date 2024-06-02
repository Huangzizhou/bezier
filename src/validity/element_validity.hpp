#pragma once
#include <vector>
#include "numbers/Interval.hpp"

namespace element_validity {
	Validity isValid(uint n, uint s, uint p, const std::vector<fp_t> &c);
}
