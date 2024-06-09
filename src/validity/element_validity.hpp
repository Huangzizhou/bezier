#pragma once
#include "transMatrices.hpp"
#include "lagrangeVector.hpp"

namespace element_validity {
	template<uint n, uint s, uint p>
	Validity isValidT(const std::vector<fp_t> &cp) {
		std::vector<Interval> out;
		lagrangeVectorT<n, s, p>(cp, out);
		for (auto x : out) std::cout << x << std::endl;
		return Validity::uncertain;
	}
}
