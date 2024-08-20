#pragma once
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include <span>
#include <array>

namespace element_validity {
	template<uint n, uint s, uint p>
	Interval lagrangeEvaluate(
		const std::span<const fp_t>,
		const std::span<const Interval>
	);
}