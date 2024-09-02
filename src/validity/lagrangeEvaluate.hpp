#pragma once
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include <array>

namespace element_validity {
	template<uint n, uint s, uint p>
	Interval lagrangeEvaluate(
		const span<const fp_t>,
		const span<const Interval>
	);
}