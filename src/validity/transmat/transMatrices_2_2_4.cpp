#include "../transMatrices.hpp"

namespace element_validity {
template<>
void initMatrices<2, 2, 4> (Matrix<Interval> &l2b, std::span<Matrix<Interval>> ssd) {
	l2b.fill({28, 0, 0, 1, 1, 1, 0, -29, 20, 1, 1, 6, 1, 1, 2, -15, 2, 1, 3, 20, 3, 1, 4, -15, 4, 1, 5, 6, 5, 1, 6, -1, 6, 2, 0, 227, 150, 2, 1, -222, 25, 2, 2, 201, 10, 2, 3, -308, 15, 2, 4, 123, 10, 2, 5, -102, 25, 2, 6, 29, 50, 3, 0, -227, 200, 3, 1, 189, 25, 3, 2, -837, 40, 3, 3, 30, 1, 3, 4, -837, 40, 3, 5, 189, 25, 3, 6, -227, 200, 4, 0, 29, 50, 4, 1, -102, 25, 4, 2, 123, 10, 4, 3, -308, 15, 4, 4, 201, 10, 4, 5, -222, 25, 4, 6, 227, 150, 5, 0, -1, 6, 5, 1, 6, 5, 5, 2, -15, 4, 5, 3, 20, 3, 5, 4, -15, 2, 5, 5, 6, 1, 5, 6, -29, 20, 6, 6, 1, 1, 7, 0, -29, 20, 7, 7, 6, 1, 7, 13, -15, 2, 7, 18, 20, 3, 7, 22, -15, 4, 7, 25, 6, 5, 7, 27, -1, 6, 8, 0, 227, 150, 8, 1, -111, 25, 8, 2, 21, 20, 8, 3, 26, 15, 8, 4, -21, 10, 8, 5, 24, 25, 8, 6, -1, 6, 8, 7, -111, 25, 8, 8, 18, 1, 8, 9, -12, 1, 8, 10, 6, 1, 8, 11, -9, 5, 8, 12, 6, 25, 8, 13, 21, 20, 8, 14, -12, 1, 8, 15, 9, 2, 8, 16, -6, 5, 8, 17, 3, 20, 8, 18, 26, 15, 8, 19, 6, 1, 8, 20, -6, 5, 8, 21, 2, 15, 8, 22, -21, 10, 8, 23, -9, 5, 8, 24, 3, 20, 8, 25, 24, 25, 8, 26, 6, 25, 8, 27, -1, 6, 9, 0, -227, 200, 9, 1, 126, 25, 9, 2, -63, 10, 9, 3, -11, 5, 9, 4, 45, 8, 9, 5, -153, 50, 9, 6, 29, 50, 9, 7, 63, 25, 9, 8, -153, 10, 9, 9, 171, 5, 9, 10, -108, 5, 9, 11, 36, 5, 9, 12, -51, 50, 9, 13, 27, 40, 9, 14, -9, 5, 9, 15, -189, 20, 9, 16, 18, 5, 9, 17, -21, 40, 9, 18, -1, 5, 9, 19, 27, 5, 9, 20, 9, 5, 9, 21, -1, 3, 9, 22, -9, 10, 9, 23, -27, 10, 9, 24, -3, 20, 9, 25, 18, 25, 9, 26, 12, 25, 9, 27, -1, 6, 10, 0, 29, 50, 10, 1, -153, 50, 10, 2, 45, 8, 10, 3, -11, 5, 10, 4, -63, 10, 10, 5, 126, 25, 10, 6, -227, 200, 10, 7, -51, 50, 10, 8, 36, 5, 10, 9, -108, 5, 10, 10, 171, 5, 10, 11, -153, 10, 10, 12, 63, 25, 10, 13, -21, 40, 10, 14, 18, 5, 10, 15, -189, 20, 10, 16, -9, 5, 10, 17, 27, 40, 10, 18, -1, 3, 10, 19, 9, 5, 10, 20, 27, 5, 10, 21, -1, 5, 10, 22, -3, 20, 10, 23, -27, 10, 10, 24, -9, 10, 10, 25, 12, 25, 10, 26, 18, 25, 10, 27, -1, 6, 11, 0, -1, 6, 11, 1, 24, 25, 11, 2, -21, 10, 11, 3, 26, 15, 11, 4, 21, 20, 11, 5, -111, 25, 11, 6, 227, 150, 11, 7, 6, 25, 11, 8, -9, 5, 11, 9, 6, 1, 11, 10, -12, 1, 11, 11, 18, 1, 11, 12, -111, 25, 11, 13, 3, 20, 11, 14, -6, 5, 11, 15, 9, 2, 11, 16, -12, 1, 11, 17, 21, 20, 11, 18, 2, 15, 11, 19, -6, 5, 11, 20, 6, 1, 11, 21, 26, 15, 11, 22, 3, 20, 11, 23, -9, 5, 11, 24, -21, 10, 11, 25, 6, 25, 11, 26, 24, 25, 11, 27, -1, 6, 12, 6, -29, 20, 12, 12, 6, 1, 12, 17, -15, 2, 12, 21, 20, 3, 12, 24, -15, 4, 12, 26, 6, 5, 12, 27, -1, 6, 13, 0, 227, 150, 13, 7, -222, 25, 13, 13, 201, 10, 13, 18, -308, 15, 13, 22, 123, 10, 13, 25, -102, 25, 13, 27, 29, 50, 14, 0, -227, 200, 14, 1, 63, 25, 14, 2, 27, 40, 14, 3, -1, 5, 14, 4, -9, 10, 14, 5, 18, 25, 14, 6, -1, 6, 14, 7, 126, 25, 14, 8, -153, 10, 14, 9, -9, 5, 14, 10, 27, 5, 14, 11, -27, 10, 14, 12, 12, 25, 14, 13, -63, 10, 14, 14, 171, 5, 14, 15, -189, 20, 14, 16, 9, 5, 14, 17, -3, 20, 14, 18, -11, 5, 14, 19, -108, 5, 14, 20, 18, 5, 14, 21, -1, 3, 14, 22, 45, 8, 14, 23, 36, 5, 14, 24, -21, 40, 14, 25, -153, 50, 14, 26, -51, 50, 14, 27, 29, 50, 15, 0, 29, 50, 15, 1, -51, 25, 15, 2, 27, 20, 15, 3, 26, 15, 15, 4, 27, 20, 15, 5, -51, 25, 15, 6, 29, 50, 15, 7, -51, 25, 15, 8, 48, 5, 15, 9, -12, 1, 15, 10, -12, 1, 15, 11, 48, 5, 15, 12, -51, 25, 15, 13, 27, 20, 15, 14, -12, 1, 15, 15, 207, 5, 15, 16, -12, 1, 15, 17, 27, 20, 15, 18, 26, 15, 15, 19, -12, 1, 15, 20, -12, 1, 15, 21, 26, 15, 15, 22, 27, 20, 15, 23, 48, 5, 15, 24, 27, 20, 15, 25, -51, 25, 15, 26, -51, 25, 15, 27, 29, 50, 16, 0, -1, 6, 16, 1, 18, 25, 16, 2, -9, 10, 16, 3, -1, 5, 16, 4, 27, 40, 16, 5, 63, 25, 16, 6, -227, 200, 16, 7, 12, 25, 16, 8, -27, 10, 16, 9, 27, 5, 16, 10, -9, 5, 16, 11, -153, 10, 16, 12, 126, 25, 16, 13, -3, 20, 16, 14, 9, 5, 16, 15, -189, 20, 16, 16, 171, 5, 16, 17, -63, 10, 16, 18, -1, 3, 16, 19, 18, 5, 16, 20, -108, 5, 16, 21, -11, 5, 16, 22, -21, 40, 16, 23, 36, 5, 16, 24, 45, 8, 16, 25, -51, 50, 16, 26, -153, 50, 16, 27, 29, 50, 17, 6, 227, 150, 17, 12, -222, 25, 17, 17, 201, 10, 17, 21, -308, 15, 17, 24, 123, 10, 17, 26, -102, 25, 17, 27, 29, 50, 18, 0, -227, 200, 18, 7, 189, 25, 18, 13, -837, 40, 18, 18, 30, 1, 18, 22, -837, 40, 18, 25, 189, 25, 18, 27, -227, 200, 19, 0, 29, 50, 19, 1, -51, 50, 19, 2, -21, 40, 19, 3, -1, 3, 19, 4, -3, 20, 19, 5, 12, 25, 19, 6, -1, 6, 19, 7, -153, 50, 19, 8, 36, 5, 19, 9, 18, 5, 19, 10, 9, 5, 19, 11, -27, 10, 19, 12, 18, 25, 19, 13, 45, 8, 19, 14, -108, 5, 19, 15, -189, 20, 19, 16, 27, 5, 19, 17, -9, 10, 19, 18, -11, 5, 19, 19, 171, 5, 19, 20, -9, 5, 19, 21, -1, 5, 19, 22, -63, 10, 19, 23, -153, 10, 19, 24, 27, 40, 19, 25, 126, 25, 19, 26, 63, 25, 19, 27, -227, 200, 20, 0, -1, 6, 20, 1, 12, 25, 20, 2, -3, 20, 20, 3, -1, 3, 20, 4, -21, 40, 20, 5, -51, 50, 20, 6, 29, 50, 20, 7, 18, 25, 20, 8, -27, 10, 20, 9, 9, 5, 20, 10, 18, 5, 20, 11, 36, 5, 20, 12, -153, 50, 20, 13, -9, 10, 20, 14, 27, 5, 20, 15, -189, 20, 20, 16, -108, 5, 20, 17, 45, 8, 20, 18, -1, 5, 20, 19, -9, 5, 20, 20, 171, 5, 20, 21, -11, 5, 20, 22, 27, 40, 20, 23, -153, 10, 20, 24, -63, 10, 20, 25, 63, 25, 20, 26, 126, 25, 20, 27, -227, 200, 21, 6, -227, 200, 21, 12, 189, 25, 21, 17, -837, 40, 21, 21, 30, 1, 21, 24, -837, 40, 21, 26, 189, 25, 21, 27, -227, 200, 22, 0, 29, 50, 22, 7, -102, 25, 22, 13, 123, 10, 22, 18, -308, 15, 22, 22, 201, 10, 22, 25, -222, 25, 22, 27, 227, 150, 23, 0, -1, 6, 23, 1, 6, 25, 23, 2, 3, 20, 23, 3, 2, 15, 23, 4, 3, 20, 23, 5, 6, 25, 23, 6, -1, 6, 23, 7, 24, 25, 23, 8, -9, 5, 23, 9, -6, 5, 23, 10, -6, 5, 23, 11, -9, 5, 23, 12, 24, 25, 23, 13, -21, 10, 23, 14, 6, 1, 23, 15, 9, 2, 23, 16, 6, 1, 23, 17, -21, 10, 23, 18, 26, 15, 23, 19, -12, 1, 23, 20, -12, 1, 23, 21, 26, 15, 23, 22, 21, 20, 23, 23, 18, 1, 23, 24, 21, 20, 23, 25, -111, 25, 23, 26, -111, 25, 23, 27, 227, 150, 24, 6, 29, 50, 24, 12, -102, 25, 24, 17, 123, 10, 24, 21, -308, 15, 24, 24, 201, 10, 24, 26, -222, 25, 24, 27, 227, 150, 25, 0, -1, 6, 25, 7, 6, 5, 25, 13, -15, 4, 25, 18, 20, 3, 25, 22, -15, 2, 25, 25, 6, 1, 25, 27, -29, 20, 26, 6, -1, 6, 26, 12, 6, 5, 26, 17, -15, 4, 26, 21, 20, 3, 26, 24, -15, 2, 26, 26, 6, 1, 26, 27, -29, 20, 27, 27, 1, 1});
	ssd[0].fill({28, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 1, 2, 2, 0, 1, 4, 2, 1, 1, 2, 2, 2, 1, 4, 3, 0, 1, 8, 3, 1, 3, 8, 3, 2, 3, 8, 3, 3, 1, 8, 4, 0, 1, 16, 4, 1, 1, 4, 4, 2, 3, 8, 4, 3, 1, 4, 4, 4, 1, 16, 5, 0, 1, 32, 5, 1, 5, 32, 5, 2, 5, 16, 5, 3, 5, 16, 5, 4, 5, 32, 5, 5, 1, 32, 6, 0, 1, 64, 6, 1, 3, 32, 6, 2, 15, 64, 6, 3, 5, 16, 6, 4, 15, 64, 6, 5, 3, 32, 6, 6, 1, 64, 7, 0, 1, 2, 7, 7, 1, 2, 8, 0, 1, 4, 8, 1, 1, 4, 8, 7, 1, 4, 8, 8, 1, 4, 9, 0, 1, 8, 9, 1, 1, 4, 9, 2, 1, 8, 9, 7, 1, 8, 9, 8, 1, 4, 9, 9, 1, 8, 10, 0, 1, 16, 10, 1, 3, 16, 10, 2, 3, 16, 10, 3, 1, 16, 10, 7, 1, 16, 10, 8, 3, 16, 10, 9, 3, 16, 10, 10, 1, 16, 11, 0, 1, 32, 11, 1, 1, 8, 11, 2, 3, 16, 11, 3, 1, 8, 11, 4, 1, 32, 11, 7, 1, 32, 11, 8, 1, 8, 11, 9, 3, 16, 11, 10, 1, 8, 11, 11, 1, 32, 12, 0, 1, 64, 12, 1, 5, 64, 12, 2, 5, 32, 12, 3, 5, 32, 12, 4, 5, 64, 12, 5, 1, 64, 12, 7, 1, 64, 12, 8, 5, 64, 12, 9, 5, 32, 12, 10, 5, 32, 12, 11, 5, 64, 12, 12, 1, 64, 13, 0, 1, 4, 13, 7, 1, 2, 13, 13, 1, 4, 14, 0, 1, 8, 14, 1, 1, 8, 14, 7, 1, 4, 14, 8, 1, 4, 14, 13, 1, 8, 14, 14, 1, 8, 15, 0, 1, 16, 15, 1, 1, 8, 15, 2, 1, 16, 15, 7, 1, 8, 15, 8, 1, 4, 15, 9, 1, 8, 15, 13, 1, 16, 15, 14, 1, 8, 15, 15, 1, 16, 16, 0, 1, 32, 16, 1, 3, 32, 16, 2, 3, 32, 16, 3, 1, 32, 16, 7, 1, 16, 16, 8, 3, 16, 16, 9, 3, 16, 16, 10, 1, 16, 16, 13, 1, 32, 16, 14, 3, 32, 16, 15, 3, 32, 16, 16, 1, 32, 17, 0, 1, 64, 17, 1, 1, 16, 17, 2, 3, 32, 17, 3, 1, 16, 17, 4, 1, 64, 17, 7, 1, 32, 17, 8, 1, 8, 17, 9, 3, 16, 17, 10, 1, 8, 17, 11, 1, 32, 17, 13, 1, 64, 17, 14, 1, 16, 17, 15, 3, 32, 17, 16, 1, 16, 17, 17, 1, 64, 18, 0, 1, 8, 18, 7, 3, 8, 18, 13, 3, 8, 18, 18, 1, 8, 19, 0, 1, 16, 19, 1, 1, 16, 19, 7, 3, 16, 19, 8, 3, 16, 19, 13, 3, 16, 19, 14, 3, 16, 19, 18, 1, 16, 19, 19, 1, 16, 20, 0, 1, 32, 20, 1, 1, 16, 20, 2, 1, 32, 20, 7, 3, 32, 20, 8, 3, 16, 20, 9, 3, 32, 20, 13, 3, 32, 20, 14, 3, 16, 20, 15, 3, 32, 20, 18, 1, 32, 20, 19, 1, 16, 20, 20, 1, 32, 21, 0, 1, 64, 21, 1, 3, 64, 21, 2, 3, 64, 21, 3, 1, 64, 21, 7, 3, 64, 21, 8, 9, 64, 21, 9, 9, 64, 21, 10, 3, 64, 21, 13, 3, 64, 21, 14, 9, 64, 21, 15, 9, 64, 21, 16, 3, 64, 21, 18, 1, 64, 21, 19, 3, 64, 21, 20, 3, 64, 21, 21, 1, 64, 22, 0, 1, 16, 22, 7, 1, 4, 22, 13, 3, 8, 22, 18, 1, 4, 22, 22, 1, 16, 23, 0, 1, 32, 23, 1, 1, 32, 23, 7, 1, 8, 23, 8, 1, 8, 23, 13, 3, 16, 23, 14, 3, 16, 23, 18, 1, 8, 23, 19, 1, 8, 23, 22, 1, 32, 23, 23, 1, 32, 24, 0, 1, 64, 24, 1, 1, 32, 24, 2, 1, 64, 24, 7, 1, 16, 24, 8, 1, 8, 24, 9, 1, 16, 24, 13, 3, 32, 24, 14, 3, 16, 24, 15, 3, 32, 24, 18, 1, 16, 24, 19, 1, 8, 24, 20, 1, 16, 24, 22, 1, 64, 24, 23, 1, 32, 24, 24, 1, 64, 25, 0, 1, 32, 25, 7, 5, 32, 25, 13, 5, 16, 25, 18, 5, 16, 25, 22, 5, 32, 25, 25, 1, 32, 26, 0, 1, 64, 26, 1, 1, 64, 26, 7, 5, 64, 26, 8, 5, 64, 26, 13, 5, 32, 26, 14, 5, 32, 26, 18, 5, 32, 26, 19, 5, 32, 26, 22, 5, 64, 26, 23, 5, 64, 26, 25, 1, 64, 26, 26, 1, 64, 27, 0, 1, 64, 27, 7, 3, 32, 27, 13, 15, 64, 27, 18, 5, 16, 27, 22, 15, 64, 27, 25, 3, 32, 27, 27, 1, 64});
	ssd[1].fill({28, 0, 0, 1, 64, 0, 7, 3, 32, 0, 13, 15, 64, 0, 18, 5, 16, 0, 22, 15, 64, 0, 25, 3, 32, 0, 27, 1, 64, 1, 1, 1, 64, 1, 7, 1, 64, 1, 8, 5, 64, 1, 13, 5, 64, 1, 14, 5, 32, 1, 18, 5, 32, 1, 19, 5, 32, 1, 22, 5, 32, 1, 23, 5, 64, 1, 25, 5, 64, 1, 26, 1, 64, 1, 27, 1, 64, 2, 2, 1, 64, 2, 8, 1, 32, 2, 9, 1, 16, 2, 13, 1, 64, 2, 14, 1, 8, 2, 15, 3, 32, 2, 18, 1, 16, 2, 19, 3, 16, 2, 20, 1, 16, 2, 22, 3, 32, 2, 23, 1, 8, 2, 24, 1, 64, 2, 25, 1, 16, 2, 26, 1, 32, 2, 27, 1, 64, 3, 3, 1, 64, 3, 9, 3, 64, 3, 10, 3, 64, 3, 14, 3, 64, 3, 15, 9, 64, 3, 16, 3, 64, 3, 18, 1, 64, 3, 19, 9, 64, 3, 20, 9, 64, 3, 21, 1, 64, 3, 22, 3, 64, 3, 23, 9, 64, 3, 24, 3, 64, 3, 25, 3, 64, 3, 26, 3, 64, 3, 27, 1, 64, 4, 4, 1, 64, 4, 10, 1, 16, 4, 11, 1, 32, 4, 15, 3, 32, 4, 16, 1, 8, 4, 17, 1, 64, 4, 19, 1, 16, 4, 20, 3, 16, 4, 21, 1, 16, 4, 22, 1, 64, 4, 23, 1, 8, 4, 24, 3, 32, 4, 25, 1, 32, 4, 26, 1, 16, 4, 27, 1, 64, 5, 5, 1, 64, 5, 11, 5, 64, 5, 12, 1, 64, 5, 16, 5, 32, 5, 17, 5, 64, 5, 20, 5, 32, 5, 21, 5, 32, 5, 23, 5, 64, 5, 24, 5, 32, 5, 25, 1, 64, 5, 26, 5, 64, 5, 27, 1, 64, 6, 6, 1, 64, 6, 12, 3, 32, 6, 17, 15, 64, 6, 21, 5, 16, 6, 24, 15, 64, 6, 26, 3, 32, 6, 27, 1, 64, 7, 7, 1, 32, 7, 13, 5, 32, 7, 18, 5, 16, 7, 22, 5, 16, 7, 25, 5, 32, 7, 27, 1, 32, 8, 8, 1, 32, 8, 13, 1, 32, 8, 14, 1, 8, 8, 18, 1, 8, 8, 19, 3, 16, 8, 22, 3, 16, 8, 23, 1, 8, 8, 25, 1, 8, 8, 26, 1, 32, 8, 27, 1, 32, 9, 9, 1, 32, 9, 14, 1, 16, 9, 15, 3, 32, 9, 18, 1, 32, 9, 19, 3, 16, 9, 20, 3, 32, 9, 22, 3, 32, 9, 23, 3, 16, 9, 24, 1, 32, 9, 25, 3, 32, 9, 26, 1, 16, 9, 27, 1, 32, 10, 10, 1, 32, 10, 15, 3, 32, 10, 16, 1, 16, 10, 19, 3, 32, 10, 20, 3, 16, 10, 21, 1, 32, 10, 22, 1, 32, 10, 23, 3, 16, 10, 24, 3, 32, 10, 25, 1, 16, 10, 26, 3, 32, 10, 27, 1, 32, 11, 11, 1, 32, 11, 16, 1, 8, 11, 17, 1, 32, 11, 20, 3, 16, 11, 21, 1, 8, 11, 23, 1, 8, 11, 24, 3, 16, 11, 25, 1, 32, 11, 26, 1, 8, 11, 27, 1, 32, 12, 12, 1, 32, 12, 17, 5, 32, 12, 21, 5, 16, 12, 24, 5, 16, 12, 26, 5, 32, 12, 27, 1, 32, 13, 13, 1, 16, 13, 18, 1, 4, 13, 22, 3, 8, 13, 25, 1, 4, 13, 27, 1, 16, 14, 14, 1, 16, 14, 18, 1, 16, 14, 19, 3, 16, 14, 22, 3, 16, 14, 23, 3, 16, 14, 25, 3, 16, 14, 26, 1, 16, 14, 27, 1, 16, 15, 15, 1, 16, 15, 19, 1, 8, 15, 20, 1, 8, 15, 22, 1, 16, 15, 23, 1, 4, 15, 24, 1, 16, 15, 25, 1, 8, 15, 26, 1, 8, 15, 27, 1, 16, 16, 16, 1, 16, 16, 20, 3, 16, 16, 21, 1, 16, 16, 23, 3, 16, 16, 24, 3, 16, 16, 25, 1, 16, 16, 26, 3, 16, 16, 27, 1, 16, 17, 17, 1, 16, 17, 21, 1, 4, 17, 24, 3, 8, 17, 26, 1, 4, 17, 27, 1, 16, 18, 18, 1, 8, 18, 22, 3, 8, 18, 25, 3, 8, 18, 27, 1, 8, 19, 19, 1, 8, 19, 22, 1, 8, 19, 23, 1, 4, 19, 25, 1, 4, 19, 26, 1, 8, 19, 27, 1, 8, 20, 20, 1, 8, 20, 23, 1, 4, 20, 24, 1, 8, 20, 25, 1, 8, 20, 26, 1, 4, 20, 27, 1, 8, 21, 21, 1, 8, 21, 24, 3, 8, 21, 26, 3, 8, 21, 27, 1, 8, 22, 22, 1, 4, 22, 25, 1, 2, 22, 27, 1, 4, 23, 23, 1, 4, 23, 25, 1, 4, 23, 26, 1, 4, 23, 27, 1, 4, 24, 24, 1, 4, 24, 26, 1, 2, 24, 27, 1, 4, 25, 25, 1, 2, 25, 27, 1, 2, 26, 26, 1, 2, 26, 27, 1, 2, 27, 27, 1, 1});
	ssd[2].fill({28, 0, 0, 1, 64, 0, 1, 3, 32, 0, 2, 15, 64, 0, 3, 5, 16, 0, 4, 15, 64, 0, 5, 3, 32, 0, 6, 1, 64, 1, 1, 1, 32, 1, 2, 5, 32, 1, 3, 5, 16, 1, 4, 5, 16, 1, 5, 5, 32, 1, 6, 1, 32, 2, 2, 1, 16, 2, 3, 1, 4, 2, 4, 3, 8, 2, 5, 1, 4, 2, 6, 1, 16, 3, 3, 1, 8, 3, 4, 3, 8, 3, 5, 3, 8, 3, 6, 1, 8, 4, 4, 1, 4, 4, 5, 1, 2, 4, 6, 1, 4, 5, 5, 1, 2, 5, 6, 1, 2, 6, 6, 1, 1, 7, 1, 1, 64, 7, 2, 5, 64, 7, 3, 5, 32, 7, 4, 5, 32, 7, 5, 5, 64, 7, 6, 1, 64, 7, 7, 1, 64, 7, 8, 5, 64, 7, 9, 5, 32, 7, 10, 5, 32, 7, 11, 5, 64, 7, 12, 1, 64, 8, 2, 1, 32, 8, 3, 1, 8, 8, 4, 3, 16, 8, 5, 1, 8, 8, 6, 1, 32, 8, 8, 1, 32, 8, 9, 1, 8, 8, 10, 3, 16, 8, 11, 1, 8, 8, 12, 1, 32, 9, 3, 1, 16, 9, 4, 3, 16, 9, 5, 3, 16, 9, 6, 1, 16, 9, 9, 1, 16, 9, 10, 3, 16, 9, 11, 3, 16, 9, 12, 1, 16, 10, 4, 1, 8, 10, 5, 1, 4, 10, 6, 1, 8, 10, 10, 1, 8, 10, 11, 1, 4, 10, 12, 1, 8, 11, 5, 1, 4, 11, 6, 1, 4, 11, 11, 1, 4, 11, 12, 1, 4, 12, 6, 1, 2, 12, 12, 1, 2, 13, 2, 1, 64, 13, 3, 1, 16, 13, 4, 3, 32, 13, 5, 1, 16, 13, 6, 1, 64, 13, 8, 1, 32, 13, 9, 1, 8, 13, 10, 3, 16, 13, 11, 1, 8, 13, 12, 1, 32, 13, 13, 1, 64, 13, 14, 1, 16, 13, 15, 3, 32, 13, 16, 1, 16, 13, 17, 1, 64, 14, 3, 1, 32, 14, 4, 3, 32, 14, 5, 3, 32, 14, 6, 1, 32, 14, 9, 1, 16, 14, 10, 3, 16, 14, 11, 3, 16, 14, 12, 1, 16, 14, 14, 1, 32, 14, 15, 3, 32, 14, 16, 3, 32, 14, 17, 1, 32, 15, 4, 1, 16, 15, 5, 1, 8, 15, 6, 1, 16, 15, 10, 1, 8, 15, 11, 1, 4, 15, 12, 1, 8, 15, 15, 1, 16, 15, 16, 1, 8, 15, 17, 1, 16, 16, 5, 1, 8, 16, 6, 1, 8, 16, 11, 1, 4, 16, 12, 1, 4, 16, 16, 1, 8, 16, 17, 1, 8, 17, 6, 1, 4, 17, 12, 1, 2, 17, 17, 1, 4, 18, 3, 1, 64, 18, 4, 3, 64, 18, 5, 3, 64, 18, 6, 1, 64, 18, 9, 3, 64, 18, 10, 9, 64, 18, 11, 9, 64, 18, 12, 3, 64, 18, 14, 3, 64, 18, 15, 9, 64, 18, 16, 9, 64, 18, 17, 3, 64, 18, 18, 1, 64, 18, 19, 3, 64, 18, 20, 3, 64, 18, 21, 1, 64, 19, 4, 1, 32, 19, 5, 1, 16, 19, 6, 1, 32, 19, 10, 3, 32, 19, 11, 3, 16, 19, 12, 3, 32, 19, 15, 3, 32, 19, 16, 3, 16, 19, 17, 3, 32, 19, 19, 1, 32, 19, 20, 1, 16, 19, 21, 1, 32, 20, 5, 1, 16, 20, 6, 1, 16, 20, 11, 3, 16, 20, 12, 3, 16, 20, 16, 3, 16, 20, 17, 3, 16, 20, 20, 1, 16, 20, 21, 1, 16, 21, 6, 1, 8, 21, 12, 3, 8, 21, 17, 3, 8, 21, 21, 1, 8, 22, 4, 1, 64, 22, 5, 1, 32, 22, 6, 1, 64, 22, 10, 1, 16, 22, 11, 1, 8, 22, 12, 1, 16, 22, 15, 3, 32, 22, 16, 3, 16, 22, 17, 3, 32, 22, 19, 1, 16, 22, 20, 1, 8, 22, 21, 1, 16, 22, 22, 1, 64, 22, 23, 1, 32, 22, 24, 1, 64, 23, 5, 1, 32, 23, 6, 1, 32, 23, 11, 1, 8, 23, 12, 1, 8, 23, 16, 3, 16, 23, 17, 3, 16, 23, 20, 1, 8, 23, 21, 1, 8, 23, 23, 1, 32, 23, 24, 1, 32, 24, 6, 1, 16, 24, 12, 1, 4, 24, 17, 3, 8, 24, 21, 1, 4, 24, 24, 1, 16, 25, 5, 1, 64, 25, 6, 1, 64, 25, 11, 5, 64, 25, 12, 5, 64, 25, 16, 5, 32, 25, 17, 5, 32, 25, 20, 5, 32, 25, 21, 5, 32, 25, 23, 5, 64, 25, 24, 5, 64, 25, 25, 1, 64, 25, 26, 1, 64, 26, 6, 1, 32, 26, 12, 5, 32, 26, 17, 5, 16, 26, 21, 5, 16, 26, 24, 5, 32, 26, 26, 1, 32, 27, 6, 1, 64, 27, 12, 3, 32, 27, 17, 15, 64, 27, 21, 5, 16, 27, 24, 15, 64, 27, 26, 3, 32, 27, 27, 1, 64});
	ssd[3].fill({28, 0, 6, 1, 64, 0, 12, 3, 32, 0, 17, 15, 64, 0, 21, 5, 16, 0, 24, 15, 64, 0, 26, 3, 32, 0, 27, 1, 64, 1, 5, 1, 64, 1, 11, 5, 64, 1, 12, 1, 64, 1, 16, 5, 32, 1, 17, 5, 64, 1, 20, 5, 32, 1, 21, 5, 32, 1, 23, 5, 64, 1, 24, 5, 32, 1, 25, 1, 64, 1, 26, 5, 64, 1, 27, 1, 64, 2, 4, 1, 64, 2, 10, 1, 16, 2, 11, 1, 32, 2, 15, 3, 32, 2, 16, 1, 8, 2, 17, 1, 64, 2, 19, 1, 16, 2, 20, 3, 16, 2, 21, 1, 16, 2, 22, 1, 64, 2, 23, 1, 8, 2, 24, 3, 32, 2, 25, 1, 32, 2, 26, 1, 16, 2, 27, 1, 64, 3, 3, 1, 64, 3, 9, 3, 64, 3, 10, 3, 64, 3, 14, 3, 64, 3, 15, 9, 64, 3, 16, 3, 64, 3, 18, 1, 64, 3, 19, 9, 64, 3, 20, 9, 64, 3, 21, 1, 64, 3, 22, 3, 64, 3, 23, 9, 64, 3, 24, 3, 64, 3, 25, 3, 64, 3, 26, 3, 64, 3, 27, 1, 64, 4, 2, 1, 64, 4, 8, 1, 32, 4, 9, 1, 16, 4, 13, 1, 64, 4, 14, 1, 8, 4, 15, 3, 32, 4, 18, 1, 16, 4, 19, 3, 16, 4, 20, 1, 16, 4, 22, 3, 32, 4, 23, 1, 8, 4, 24, 1, 64, 4, 25, 1, 16, 4, 26, 1, 32, 4, 27, 1, 64, 5, 1, 1, 64, 5, 7, 1, 64, 5, 8, 5, 64, 5, 13, 5, 64, 5, 14, 5, 32, 5, 18, 5, 32, 5, 19, 5, 32, 5, 22, 5, 32, 5, 23, 5, 64, 5, 25, 5, 64, 5, 26, 1, 64, 5, 27, 1, 64, 6, 0, 1, 64, 6, 7, 3, 32, 6, 13, 15, 64, 6, 18, 5, 16, 6, 22, 15, 64, 6, 25, 3, 32, 6, 27, 1, 64, 7, 5, 1, 64, 7, 6, 1, 64, 7, 11, 5, 64, 7, 12, 5, 64, 7, 16, 5, 32, 7, 17, 5, 32, 7, 20, 5, 32, 7, 21, 5, 32, 7, 23, 5, 64, 7, 24, 5, 64, 7, 25, 1, 64, 7, 26, 1, 64, 8, 4, 1, 64, 8, 5, 1, 64, 8, 10, 1, 16, 8, 11, 5, 64, 8, 12, 1, 64, 8, 15, 3, 32, 8, 16, 5, 32, 8, 17, 1, 16, 8, 19, 1, 16, 8, 20, 5, 32, 8, 21, 3, 32, 8, 22, 1, 64, 8, 23, 5, 64, 8, 24, 1, 16, 8, 25, 1, 64, 8, 26, 1, 64, 9, 3, 1, 64, 9, 4, 1, 64, 9, 9, 3, 64, 9, 10, 5, 64, 9, 11, 1, 32, 9, 14, 3, 64, 9, 15, 9, 64, 9, 16, 7, 64, 9, 17, 1, 64, 9, 18, 1, 64, 9, 19, 7, 64, 9, 20, 9, 64, 9, 21, 3, 64, 9, 22, 1, 32, 9, 23, 5, 64, 9, 24, 3, 64, 9, 25, 1, 64, 9, 26, 1, 64, 10, 2, 1, 64, 10, 3, 1, 64, 10, 8, 1, 32, 10, 9, 5, 64, 10, 10, 3, 64, 10, 13, 1, 64, 10, 14, 7, 64, 10, 15, 9, 64, 10, 16, 3, 64, 10, 18, 3, 64, 10, 19, 9, 64, 10, 20, 7, 64, 10, 21, 1, 64, 10, 22, 3, 64, 10, 23, 5, 64, 10, 24, 1, 32, 10, 25, 1, 64, 10, 26, 1, 64, 11, 1, 1, 64, 11, 2, 1, 64, 11, 7, 1, 64, 11, 8, 5, 64, 11, 9, 1, 16, 11, 13, 1, 16, 11, 14, 5, 32, 11, 15, 3, 32, 11, 18, 3, 32, 11, 19, 5, 32, 11, 20, 1, 16, 11, 22, 1, 16, 11, 23, 5, 64, 11, 24, 1, 64, 11, 25, 1, 64, 11, 26, 1, 64, 12, 0, 1, 64, 12, 1, 1, 64, 12, 7, 5, 64, 12, 8, 5, 64, 12, 13, 5, 32, 12, 14, 5, 32, 12, 18, 5, 32, 12, 19, 5, 32, 12, 22, 5, 64, 12, 23, 5, 64, 12, 25, 1, 64, 12, 26, 1, 64, 13, 4, 1, 64, 13, 5, 1, 32, 13, 6, 1, 64, 13, 10, 1, 16, 13, 11, 1, 8, 13, 12, 1, 16, 13, 15, 3, 32, 13, 16, 3, 16, 13, 17, 3, 32, 13, 19, 1, 16, 13, 20, 1, 8, 13, 21, 1, 16, 13, 22, 1, 64, 13, 23, 1, 32, 13, 24, 1, 64, 14, 3, 1, 64, 14, 4, 1, 32, 14, 5, 1, 64, 14, 9, 3, 64, 14, 10, 7, 64, 14, 11, 5, 64, 14, 12, 1, 64, 14, 14, 3, 64, 14, 15, 9, 64, 14, 16, 9, 64, 14, 17, 3, 64, 14, 18, 1, 64, 14, 19, 5, 64, 14, 20, 7, 64, 14, 21, 3, 64, 14, 22, 1, 64, 14, 23, 1, 32, 14, 24, 1, 64, 15, 2, 1, 64, 15, 3, 1, 32, 15, 4, 1, 64, 15, 8, 1, 32, 15, 9, 3, 32, 15, 10, 3, 32, 15, 11, 1, 32, 15, 13, 1, 64, 15, 14, 3, 32, 15, 15, 5, 32, 15, 16, 3, 32, 15, 17, 1, 64, 15, 18, 1, 32, 15, 19, 3, 32, 15, 20, 3, 32, 15, 21, 1, 32, 15, 22, 1, 64, 15, 23, 1, 32, 15, 24, 1, 64, 16, 1, 1, 64, 16, 2, 1, 32, 16, 3, 1, 64, 16, 7, 1, 64, 16, 8, 5, 64, 16, 9, 7, 64, 16, 10, 3, 64, 16, 13, 3, 64, 16, 14, 9, 64, 16, 15, 9, 64, 16, 16, 3, 64, 16, 18, 3, 64, 16, 19, 7, 64, 16, 20, 5, 64, 16, 21, 1, 64, 16, 22, 1, 64, 16, 23, 1, 32, 16, 24, 1, 64, 17, 0, 1, 64, 17, 1, 1, 32, 17, 2, 1, 64, 17, 7, 1, 16, 17, 8, 1, 8, 17, 9, 1, 16, 17, 13, 3, 32, 17, 14, 3, 16, 17, 15, 3, 32, 17, 18, 1, 16, 17, 19, 1, 8, 17, 20, 1, 16, 17, 22, 1, 64, 17, 23, 1, 32, 17, 24, 1, 64, 18, 3, 1, 64, 18, 4, 3, 64, 18, 5, 3, 64, 18, 6, 1, 64, 18, 9, 3, 64, 18, 10, 9, 64, 18, 11, 9, 64, 18, 12, 3, 64, 18, 14, 3, 64, 18, 15, 9, 64, 18, 16, 9, 64, 18, 17, 3, 64, 18, 18, 1, 64, 18, 19, 3, 64, 18, 20, 3, 64, 18, 21, 1, 64, 19, 2, 1, 64, 19, 3, 3, 64, 19, 4, 3, 64, 19, 5, 1, 64, 19, 8, 1, 32, 19, 9, 7, 64, 19, 10, 9, 64, 19, 11, 5, 64, 19, 12, 1, 64, 19, 13, 1, 64, 19, 14, 5, 64, 19, 15, 9, 64, 19, 16, 7, 64, 19, 17, 1, 32, 19, 18, 1, 64, 19, 19, 3, 64, 19, 20, 3, 64, 19, 21, 1, 64, 20, 1, 1, 64, 20, 2, 3, 64, 20, 3, 3, 64, 20, 4, 1, 64, 20, 7, 1, 64, 20, 8, 5, 64, 20, 9, 9, 64, 20, 10, 7, 64, 20, 11, 1, 32, 20, 13, 1, 32, 20, 14, 7, 64, 20, 15, 9, 64, 20, 16, 5, 64, 20, 17, 1, 64, 20, 18, 1, 64, 20, 19, 3, 64, 20, 20, 3, 64, 20, 21, 1, 64, 21, 0, 1, 64, 21, 1, 3, 64, 21, 2, 3, 64, 21, 3, 1, 64, 21, 7, 3, 64, 21, 8, 9, 64, 21, 9, 9, 64, 21, 10, 3, 64, 21, 13, 3, 64, 21, 14, 9, 64, 21, 15, 9, 64, 21, 16, 3, 64, 21, 18, 1, 64, 21, 19, 3, 64, 21, 20, 3, 64, 21, 21, 1, 64, 22, 2, 1, 64, 22, 3, 1, 16, 22, 4, 3, 32, 22, 5, 1, 16, 22, 6, 1, 64, 22, 8, 1, 32, 22, 9, 1, 8, 22, 10, 3, 16, 22, 11, 1, 8, 22, 12, 1, 32, 22, 13, 1, 64, 22, 14, 1, 16, 22, 15, 3, 32, 22, 16, 1, 16, 22, 17, 1, 64, 23, 1, 1, 64, 23, 2, 1, 16, 23, 3, 3, 32, 23, 4, 1, 16, 23, 5, 1, 64, 23, 7, 1, 64, 23, 8, 5, 64, 23, 9, 5, 32, 23, 10, 5, 32, 23, 11, 5, 64, 23, 12, 1, 64, 23, 13, 1, 64, 23, 14, 1, 16, 23, 15, 3, 32, 23, 16, 1, 16, 23, 17, 1, 64, 24, 0, 1, 64, 24, 1, 1, 16, 24, 2, 3, 32, 24, 3, 1, 16, 24, 4, 1, 64, 24, 7, 1, 32, 24, 8, 1, 8, 24, 9, 3, 16, 24, 10, 1, 8, 24, 11, 1, 32, 24, 13, 1, 64, 24, 14, 1, 16, 24, 15, 3, 32, 24, 16, 1, 16, 24, 17, 1, 64, 25, 1, 1, 64, 25, 2, 5, 64, 25, 3, 5, 32, 25, 4, 5, 32, 25, 5, 5, 64, 25, 6, 1, 64, 25, 7, 1, 64, 25, 8, 5, 64, 25, 9, 5, 32, 25, 10, 5, 32, 25, 11, 5, 64, 25, 12, 1, 64, 26, 0, 1, 64, 26, 1, 5, 64, 26, 2, 5, 32, 26, 3, 5, 32, 26, 4, 5, 64, 26, 5, 1, 64, 26, 7, 1, 64, 26, 8, 5, 64, 26, 9, 5, 32, 26, 10, 5, 32, 26, 11, 5, 64, 26, 12, 1, 64, 27, 0, 1, 64, 27, 1, 3, 32, 27, 2, 15, 64, 27, 3, 5, 16, 27, 4, 15, 64, 27, 5, 3, 32, 27, 6, 1, 64});
}
}