#include "transMatrices.hpp"

namespace element_validity {
template<>
void initMatricesT<1, 1, 5> (Matrix<Interval> &l2b, std::pair<Matrix<Interval>, Matrix<Interval>> &tsd, std::array<Matrix<Interval>, 4> &ssd) {
	l2b.fill({10, 0, 0, 1, 1, 1, 1, 1, 1, 2, 0, -13, 12, 2, 2, 4, 1, 2, 4, -3, 1, 2, 6, 4, 3, 2, 8, -1, 4, 3, 1, -13, 12, 3, 3, 4, 1, 3, 5, -3, 1, 3, 7, 4, 3, 3, 9, -1, 4, 4, 0, 13, 18, 4, 2, -32, 9, 4, 4, 20, 3, 4, 6, -32, 9, 4, 8, 13, 18, 5, 1, 13, 18, 5, 3, -32, 9, 5, 5, 20, 3, 5, 7, -32, 9, 5, 9, 13, 18, 6, 0, -1, 4, 6, 2, 4, 3, 6, 4, -3, 1, 6, 6, 4, 1, 6, 8, -13, 12, 7, 1, -1, 4, 7, 3, 4, 3, 7, 5, -3, 1, 7, 7, 4, 1, 7, 9, -13, 12, 8, 8, 1, 1, 9, 9, 1, 1});
	tsd.first.fill({10, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2, 2, 2, 1, 1, 3, 2, 1, 2, 3, 3, 1, 2, 4, 4, 1, 1, 5, 4, 1, 2, 5, 5, 1, 2, 6, 6, 1, 1, 7, 6, 1, 2, 7, 7, 1, 2, 8, 8, 1, 1, 9, 8, 1, 2, 9, 9, 1, 2});
	tsd.second.fill({10, 0, 0, 1, 2, 0, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2, 2, 3, 1, 2, 3, 3, 1, 1, 4, 4, 1, 2, 4, 5, 1, 2, 5, 5, 1, 1, 6, 6, 1, 2, 6, 7, 1, 2, 7, 7, 1, 1, 8, 8, 1, 2, 8, 9, 1, 2, 9, 9, 1, 1});
	ssd[0].fill({10, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2, 2, 0, 1, 2, 2, 2, 1, 2, 3, 0, 1, 4, 3, 1, 1, 4, 3, 2, 1, 4, 3, 3, 1, 4, 4, 0, 1, 4, 4, 2, 1, 2, 4, 4, 1, 4, 5, 0, 1, 8, 5, 1, 1, 8, 5, 2, 1, 4, 5, 3, 1, 4, 5, 4, 1, 8, 5, 5, 1, 8, 6, 0, 1, 8, 6, 2, 3, 8, 6, 4, 3, 8, 6, 6, 1, 8, 7, 0, 1, 16, 7, 1, 1, 16, 7, 2, 3, 16, 7, 3, 3, 16, 7, 4, 3, 16, 7, 5, 3, 16, 7, 6, 1, 16, 7, 7, 1, 16, 8, 0, 1, 16, 8, 2, 1, 4, 8, 4, 3, 8, 8, 6, 1, 4, 8, 8, 1, 16, 9, 0, 1, 32, 9, 1, 1, 32, 9, 2, 1, 8, 9, 3, 1, 8, 9, 4, 3, 16, 9, 5, 3, 16, 9, 6, 1, 8, 9, 7, 1, 8, 9, 8, 1, 32, 9, 9, 1, 32});
	ssd[1].fill({10, 0, 0, 1, 16, 0, 2, 1, 4, 0, 4, 3, 8, 0, 6, 1, 4, 0, 8, 1, 16, 1, 0, 1, 32, 1, 1, 1, 32, 1, 2, 1, 8, 1, 3, 1, 8, 1, 4, 3, 16, 1, 5, 3, 16, 1, 6, 1, 8, 1, 7, 1, 8, 1, 8, 1, 32, 1, 9, 1, 32, 2, 2, 1, 8, 2, 4, 3, 8, 2, 6, 3, 8, 2, 8, 1, 8, 3, 2, 1, 16, 3, 3, 1, 16, 3, 4, 3, 16, 3, 5, 3, 16, 3, 6, 3, 16, 3, 7, 3, 16, 3, 8, 1, 16, 3, 9, 1, 16, 4, 4, 1, 4, 4, 6, 1, 2, 4, 8, 1, 4, 5, 4, 1, 8, 5, 5, 1, 8, 5, 6, 1, 4, 5, 7, 1, 4, 5, 8, 1, 8, 5, 9, 1, 8, 6, 6, 1, 2, 6, 8, 1, 2, 7, 6, 1, 4, 7, 7, 1, 4, 7, 8, 1, 4, 7, 9, 1, 4, 8, 8, 1, 1, 9, 8, 1, 2, 9, 9, 1, 2});
	ssd[2].fill({10, 0, 0, 1, 2, 0, 1, 1, 2, 1, 1, 1, 1, 2, 0, 1, 4, 2, 1, 1, 4, 2, 2, 1, 4, 2, 3, 1, 4, 3, 1, 1, 2, 3, 3, 1, 2, 4, 0, 1, 8, 4, 1, 1, 8, 4, 2, 1, 4, 4, 3, 1, 4, 4, 4, 1, 8, 4, 5, 1, 8, 5, 1, 1, 4, 5, 3, 1, 2, 5, 5, 1, 4, 6, 0, 1, 16, 6, 1, 1, 16, 6, 2, 3, 16, 6, 3, 3, 16, 6, 4, 3, 16, 6, 5, 3, 16, 6, 6, 1, 16, 6, 7, 1, 16, 7, 1, 1, 8, 7, 3, 3, 8, 7, 5, 3, 8, 7, 7, 1, 8, 8, 0, 1, 32, 8, 1, 1, 32, 8, 2, 1, 8, 8, 3, 1, 8, 8, 4, 3, 16, 8, 5, 3, 16, 8, 6, 1, 8, 8, 7, 1, 8, 8, 8, 1, 32, 8, 9, 1, 32, 9, 1, 1, 16, 9, 3, 1, 4, 9, 5, 3, 8, 9, 7, 1, 4, 9, 9, 1, 16});
	ssd[3].fill({10, 0, 0, 1, 32, 0, 1, 1, 32, 0, 2, 1, 8, 0, 3, 1, 8, 0, 4, 3, 16, 0, 5, 3, 16, 0, 6, 1, 8, 0, 7, 1, 8, 0, 8, 1, 32, 0, 9, 1, 32, 1, 1, 1, 16, 1, 3, 1, 4, 1, 5, 3, 8, 1, 7, 1, 4, 1, 9, 1, 16, 2, 2, 1, 16, 2, 3, 1, 16, 2, 4, 3, 16, 2, 5, 3, 16, 2, 6, 3, 16, 2, 7, 3, 16, 2, 8, 1, 16, 2, 9, 1, 16, 3, 3, 1, 8, 3, 5, 3, 8, 3, 7, 3, 8, 3, 9, 1, 8, 4, 4, 1, 8, 4, 5, 1, 8, 4, 6, 1, 4, 4, 7, 1, 4, 4, 8, 1, 8, 4, 9, 1, 8, 5, 5, 1, 4, 5, 7, 1, 2, 5, 9, 1, 4, 6, 6, 1, 4, 6, 7, 1, 4, 6, 8, 1, 4, 6, 9, 1, 4, 7, 7, 1, 2, 7, 9, 1, 2, 8, 8, 1, 2, 8, 9, 1, 2, 9, 9, 1, 1});
}
}