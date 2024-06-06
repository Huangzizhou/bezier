#pragma once
#include <array>
#include "numbers/Interval.hpp"
#include "numbers/Rational.hpp"

namespace element_validity {
	template<uint n, uint s, uint p, bool t>
	void lagrangeVector(const std::vector<fp_t>&, std::vector<Interval>&);

	template<uint n, uint s, uint p, bool t>
	Validity isValid(const std::vector<fp_t> &cp) {
		std::vector<Interval> out;
		lagrangeVector<n, s, p, t>(cp, out);
		for (auto x : out) std::cout << x << std::endl;
		return Validity::uncertain;
	}
}
