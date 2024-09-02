#pragma once
#include "numbers/Matrix.hpp"

namespace element_validity {
	template<uint n, uint s, uint p>
	void initMatrices(
		Matrix<Interval> &l2b,
		span<Matrix<Interval>> ssd
	);
	
	template<uint n, uint s, uint p>
	void initMatricesT(
		Matrix<Interval> &l2b,
		span<Matrix<Interval>> tsd,
		span<Matrix<Interval>> ssd
	);
}