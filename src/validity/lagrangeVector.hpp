#pragma once
#include "numbers/Interval.hpp"
#include <array>

namespace element_validity {
	template<uint n, uint s, uint p>
	void lagrangeVector(const span<const fp_t>, const span<Interval>);

	template<uint n, uint s, uint p>
	void lagrangeVectorT(const span<const fp_t>, const span<Interval>);
}