#pragma once
#include "numbers/Interval.hpp"

namespace element_validity {
	template<uint n, uint s, uint p>
	void lagrangeVectorT(const std::vector<fp_t>&, std::vector<Interval>&);
}