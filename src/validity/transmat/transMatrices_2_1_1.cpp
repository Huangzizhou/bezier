#include "../transMatrices.hpp"

namespace element_validity {
template<>
void initMatrices<2, 1, 1> (Matrix<Interval> &l2b, span<Matrix<Interval>> ssd) {
	l2b.fill({4, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 3, 3, 1, 1});
	ssd[0].fill({4, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2, 2, 0, 1, 2, 2, 2, 1, 2, 3, 0, 1, 4, 3, 1, 1, 4, 3, 2, 1, 4, 3, 3, 1, 4});
	ssd[1].fill({4, 0, 0, 1, 2, 0, 2, 1, 2, 1, 0, 1, 4, 1, 1, 1, 4, 1, 2, 1, 4, 1, 3, 1, 4, 2, 2, 1, 1, 3, 2, 1, 2, 3, 3, 1, 2});
	ssd[2].fill({4, 0, 0, 1, 2, 0, 1, 1, 2, 1, 1, 1, 1, 2, 0, 1, 4, 2, 1, 1, 4, 2, 2, 1, 4, 2, 3, 1, 4, 3, 1, 1, 2, 3, 3, 1, 2});
	ssd[3].fill({4, 0, 0, 1, 4, 0, 1, 1, 4, 0, 2, 1, 4, 0, 3, 1, 4, 1, 1, 1, 2, 1, 3, 1, 2, 2, 2, 1, 2, 2, 3, 1, 2, 3, 3, 1, 1});
}
}
