#include "transMatrices.hpp"

namespace element_validity {
template<>
void initMatricesT<1, 1, 1> (Matrix<Interval> &l2b, std::span<Matrix<Interval>> tsd, std::span<Matrix<Interval>> ssd) {
	l2b.fill({2, 0, 0, 1, 1, 1, 1, 1, 1});
	tsd[0].fill({2, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2});
	tsd[1].fill({2, 0, 0, 1, 2, 0, 1, 1, 2, 1, 1, 1, 1});
	ssd[0].fill({2, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2});
	ssd[1].fill({2, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2});
	ssd[2].fill({2, 0, 0, 1, 2, 0, 1, 1, 2, 1, 1, 1, 1});
	ssd[3].fill({2, 0, 0, 1, 2, 0, 1, 1, 2, 1, 1, 1, 1});
}
}
