#pragma once
#include "numbers/Matrix.hpp"
#include <span>

namespace element_validity {
	template<uint n, uint s, uint p>
	void initMatrices(
		Matrix<Interval> &l2b,
		std::span<Matrix<Interval>> ssd
	);
	
	template<uint n, uint s, uint p>
	void initMatricesT(
		Matrix<Interval> &l2b,
		std::span<Matrix<Interval>> tsd,
		std::span<Matrix<Interval>> ssd
	);
}