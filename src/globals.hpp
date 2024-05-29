#pragma once
#include <limits>
#include <iomanip>
#include <cassert>

/*
Things that should be available everywhere in the project
*/

namespace element_validity {
	// unsigned integer
	using uint = unsigned;

	// floating point type
	using fp_t = double;
	const auto fp_fmt =
		std::setprecision(std::numeric_limits<fp_t>::max_digits10);

	// type for timers
	using chrono_t = double;
}