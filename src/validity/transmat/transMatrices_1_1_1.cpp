#include "../transMatrices.hpp"

namespace element_validity {
template<>
void initMatrices<1, 1, 1> (Matrix<Interval> &l2b, std::span<Matrix<Interval>> ssd) {
	l2b.fill({1, 0, 0, 1, 1});
	ssd[0].fill({1, 0, 0, 1, 1});
	ssd[1].fill({1, 0, 0, 1, 1});
}
}
