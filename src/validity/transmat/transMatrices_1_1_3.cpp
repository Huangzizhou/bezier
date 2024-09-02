#include "../transMatrices.hpp"

namespace element_validity {
template<>
void initMatrices<1, 1, 3> (Matrix<Interval> &l2b, span<Matrix<Interval>> ssd) {
	l2b.fill({3, 0, 0, 1, 1, 1, 0, -1, 2, 1, 1, 2, 1, 1, 2, -1, 2, 2, 2, 1, 1});
	ssd[0].fill({3, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2, 2, 0, 1, 4, 2, 1, 1, 2, 2, 2, 1, 4});
	ssd[1].fill({3, 0, 0, 1, 4, 0, 1, 1, 2, 0, 2, 1, 4, 1, 1, 1, 2, 1, 2, 1, 2, 2, 2, 1, 1});
}
}
