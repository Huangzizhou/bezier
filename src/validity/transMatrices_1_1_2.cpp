#include "transMatrices.hpp"

namespace element_validity {
template<>
void initMatrices<1, 1, 2> (Matrix<Interval> &l2b, std::span<Matrix<Interval>> ssd) {
	l2b.fill({2, 0, 0, 1, 1, 1, 1, 1, 1});
	ssd[0].fill({2, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2});
	ssd[1].fill({2, 0, 0, 1, 2, 0, 1, 1, 2, 1, 1, 1, 1});
}
}
