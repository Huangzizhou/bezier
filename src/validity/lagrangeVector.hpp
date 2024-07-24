#pragma once
#include "numbers/Interval.hpp"
#include <span>

namespace element_validity {
	template<uint n, uint s, uint p>
	void lagrangeVector(const std::span<const fp_t>, const std::span<Interval>);

	template<uint n, uint s, uint p>
	void lagrangeVectorT(const std::span<const fp_t>, const std::span<Interval>);
}