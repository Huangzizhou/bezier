#include "../lagrangeEvaluate.hpp"

#define R(p, q) (Interval(p) / q)
#define I const Interval 

namespace element_validity {
template<>
Interval lagrangeEvaluate<1, 1, 5>(
	const std::span<const fp_t> xFP,
	const std::span<const Interval> lagVec
) {
	assert(xFP.size() == 1);
	assert(lagVec.size() == 6);
	std::array<Interval, 1> x;
	for (uint i = 0; i < 1; ++i) x[i] = xFP[i];
	Interval acc = 0.;
	I tmp_0 = powi(x[0], 2);
	I tmp_1 = powi(x[0], 3);
	I tmp_2 = powi(x[0], 4);
	I tmp_3 = powi(x[0], 5);
	I tmp_4 = (R(625, 24))*tmp_3;
	I tmp_5 = 25*x[0];
	I tmp_6 = (R(3125, 24))*tmp_3;
	I tmp_7 = (R(3125, 12))*tmp_3;
	acc += lagVec[0] * (R(375, 8))*tmp_0 - R(2125, 24)*tmp_1 + (R(625, 8))*tmp_2 - tmp_4 - R(137, 12)*x[0] + 1;
	acc += lagVec[1] * -R(1925, 12)*tmp_0 + (R(8875, 24))*tmp_1 - R(4375, 12)*tmp_2 + tmp_5 + tmp_6;
	acc += lagVec[2] * (R(2675, 12))*tmp_0 - R(7375, 12)*tmp_1 + (R(8125, 12))*tmp_2 - tmp_5 - tmp_7;
	acc += lagVec[3] * -R(325, 2)*tmp_0 + (R(6125, 12))*tmp_1 - 625*tmp_2 + tmp_7 + (R(50, 3))*x[0];
	acc += lagVec[4] * (R(1525, 24))*tmp_0 - R(5125, 24)*tmp_1 + (R(6875, 24))*tmp_2 - tmp_6 - R(25, 4)*x[0];
	acc += lagVec[5] * -R(125, 12)*tmp_0 + (R(875, 24))*tmp_1 - R(625, 12)*tmp_2 + tmp_4 + x[0];
	return acc;
}}
#undef R
#undef I