#include "../transMatrices.hpp"

namespace element_validity {
template<>
void initMatrices<2, 2, 1> (Matrix<Interval> &l2b, span<Matrix<Interval>> ssd) {
	l2b.fill({1, 0, 0, 1, 1});
	ssd[0].fill({1, 0, 0, 1, 1});
	ssd[1].fill({1, 0, 0, 1, 1});
	ssd[2].fill({1, 0, 0, 1, 1});
	ssd[3].fill({1, 0, 0, 1, 1});
}
}
