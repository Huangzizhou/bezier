#pragma once
#include "numbers/Matrix.hpp"
#include <array>

namespace element_validity {
	template<uint n, uint s, uint p>
	void initMatricesT(
		Matrix<Interval> &l2b,
		std::pair<Matrix<Interval>, Matrix<Interval>> &tsd,
		std::array<Matrix<Interval>, 2<<n> &ssd
	);
}