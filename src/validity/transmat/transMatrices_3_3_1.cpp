#include "../transMatrices.hpp"

namespace element_validity {
template<>
void initMatrices<3, 3, 1> (Matrix<Interval> &l2b, span<Matrix<Interval>> ssd) {
	l2b.fill({1, 0, 0, 1, 1});
	ssd[0].fill({1, 0, 0, 1, 1});
	ssd[1].fill({1, 0, 0, 1, 1});
	ssd[2].fill({1, 0, 0, 1, 1});
	ssd[3].fill({1, 0, 0, 1, 1});
	ssd[4].fill({1, 0, 0, 1, 1});
	ssd[5].fill({1, 0, 0, 1, 1});
	ssd[6].fill({1, 0, 0, 1, 1});
	ssd[7].fill({1, 0, 0, 1, 1});
}
}
