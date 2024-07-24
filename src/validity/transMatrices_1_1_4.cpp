#include "transMatrices.hpp"

namespace element_validity {
template<>
void initMatrices<1, 1, 4> (Matrix<Interval> &l2b, std::span<Matrix<Interval>> ssd) {
	l2b.fill({4, 0, 0, 1, 1, 1, 0, -5, 6, 1, 1, 3, 1, 1, 2, -3, 2, 1, 3, 1, 3, 2, 0, 1, 3, 2, 1, -3, 2, 2, 2, 3, 1, 2, 3, -5, 6, 3, 3, 1, 1});
	ssd[0].fill({4, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2, 2, 0, 1, 4, 2, 1, 1, 2, 2, 2, 1, 4, 3, 0, 1, 8, 3, 1, 3, 8, 3, 2, 3, 8, 3, 3, 1, 8});
	ssd[1].fill({4, 0, 0, 1, 8, 0, 1, 3, 8, 0, 2, 3, 8, 0, 3, 1, 8, 1, 1, 1, 4, 1, 2, 1, 2, 1, 3, 1, 4, 2, 2, 1, 2, 2, 3, 1, 2, 3, 3, 1, 1});
}
}
